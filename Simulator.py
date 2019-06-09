import pysam, random, time, sys, numpy, csv


def is_a_valid_read(chr):
    '''
    check if the chromosome name is valid or not
    :param chr:
    :return:
    '''
    if len(chr) < 3 or chr[0:3] != 'chr':
        return False

    if len(chr) == 4 and (chr[3].upper() == 'X' or chr[3].upper() == 'Y'):
        return True
    try:
        if int(chr[3:]) > 0 and int(chr[3:]) < 24:
            return True
    except Exception:
        return False


def has_intersection(begin1, end1, begin2, end2):
    '''
    check if two ranges have intersection with each other or not
    :param begin1:
    :param end1:
    :param begin2:
    :param end2:
    :return:
    '''
    try:
        if begin1 <= begin2 and end1 >= begin2:
            return True
        if begin2 <= begin1 and end2 >= begin1:
            return True
        return False
    except Exception as e:
        print(e)
        print(begin1, " ", end1, " ", begin2, " ", end2)
        exit(1)


def get_read_in_range(all_reads, window_size, list_of_n, p_tumor, list_of_positions_errors):
    '''

    :param all_reads: list of [ type (tumor or normal ), chromosome, begin_pos, end_pos, raw_read_string]
    :param begin: first position in this chromosome
    :param window_size:
    :param list_of_n:
    :param p_tumor
    :param p_normal
    :param list_of_positions_errors: list of [pos, normal error rate, tissue error rate]
    :return:
    '''

    reads = []
    normal_simulated = 0
    tissue_simulated = 0
    new_position_error = []
    if list_of_n != None:
        while len(all_reads) > 0:
            p_tumor_val = p_tumor.pop(0)
            p_normal = 1 - p_tumor_val
            i = int(all_reads[0][2] / window_size)
            n = list_of_n.pop(0)
            begin_pos = i * window_size
            end_pos = begin_pos + window_size
            window_read_tumor = []
            window_read_normal = []
            if len(list_of_positions_errors) > 0 and begin_pos <= list_of_positions_errors[0][0] <= end_pos:
                list_val = list_of_positions_errors.pop(0)
                new_position_error.append([list_val[0], p_normal * list_val[1] + p_tumor_val * list_val[2]])
            while len(all_reads) > 0 and all_reads[0][3] < begin_pos:
                all_reads.pop(0)
            while len(all_reads) > 0 and all_reads[0][2] <= end_pos:
                read = all_reads.pop(0)
                if read[0] == 'tumor':
                    window_read_tumor.append(read[4])
                elif read[0] == 'normal':
                    window_read_normal.append(read[4])
            # n = len(window_read_tumor) + len(window_read_normal)
            if len(window_read_normal) == 0:
                selected_normal_reads = []
            else:
                selected_normal_reads = numpy.random.choice(window_read_normal,
                                                            min(int(n * p_normal), len(window_read_normal)),
                                                            replace=False)
            if len(window_read_tumor) == 0:
                selected_tumor_reads = []
            else:
                selected_tumor_reads = numpy.random.choice(window_read_tumor,
                                                           min(int(n * p_tumor_val), len(window_read_tumor)),
                                                           replace=False)
            for r in selected_normal_reads:
                reads.append(r)
            for r in selected_tumor_reads:
                reads.append(r)
            normal_simulated += len(selected_normal_reads)
            tissue_simulated += len(selected_tumor_reads)
        return add_PCR_error(reads, new_position_error), normal_simulated, tissue_simulated

    # else:
    # p_normal = 1 - p_tumor
    # n = 300
    # while len(all_reads) > 0:
    #     i = int(all_reads[0][2] / window_size)
    #     begin_pos =  i * window_size
    #     end_pos = begin_pos + window_size
    #     window_read_tumor = []
    #     window_read_normal = []
    #     while len(all_reads) > 0 and all_reads[0][3] < begin_pos:
    #         all_reads.pop(0)
    #     while len(all_reads) > 0 and all_reads[0][2] <= end_pos:
    #         read = all_reads.pop(0)
    #         if read[0] == 'tumor':
    #             window_read_tumor.append(read[4])
    #         elif read[0] == 'normal' and random.random() < p_normal:
    #             window_read_normal.append(read[4])
    #     if len(window_read_normal) == 0:
    #         selected_normal_reads = []
    #     else:
    #         selected_normal_reads = numpy.random.choice(window_read_normal,
    #                                                     min(int(n * p_normal), len(window_read_normal)),
    #                                                     replace=False)
    #     if len(window_read_tumor) == 0:
    #         selected_tumor_reads = []
    #     else:
    #         selected_tumor_reads = numpy.random.choice(window_read_tumor,
    #                                                    min(int(n * p_tumor), len(window_read_tumor)),
    #                                                    replace=False)
    #     for r in selected_normal_reads:
    #         reads.append(r)
    #     for r in selected_tumor_reads:
    #         reads.append(r)

    return reads


def is_in_matched_region(position, cigartuples):
    matched = 0
    end = cigartuples[0][1]
    while position > end and len(cigartuples) > 0:
        end = end + cigartuples[0][1]
        if cigartuples[0][0] == 0:
            matched += cigartuples[0][1]
        cigartuples.pop(0)
    if len(cigartuples) == 0:
        return False, None
    if cigartuples[0][0] == 0:
        return True, matched
    return False, None


def add_PCR_error(all_reads, list_of_positions_errors):
    '''
        :param all_reads:
        :param list_of_positions: list of sorted all positions we want to add pcr error
        :return:
        '''
    reads = []
    for read in all_reads:
        read_begin = read.get_reference_positions()[0]
        read_end = read.get_reference_positions()[-1]
        read_object = [read, False, []]
        for val in list_of_positions_errors:
            pos = val[0]
            p_error = val[1]
            if pos < read_begin:
                continue
            if pos > read_end:
                break
            if random.random() < p_error:
                rel_pos = pos - read_begin
                is_in_region, matched = is_in_matched_region(rel_pos, read.cigartuples.copy())
                if is_in_region:
                    nucs = ['A', 'C', 'T', 'G']
                    new_nuc = numpy.random.choice(nucs, 1)[0]
                    read_object[1] = True
                    read_object[2].append([rel_pos - matched, new_nuc])
        reads.append(read_object)
    return reads


def read_to_string(read, rNext, pNext):
    '''
    convert pysam read object to samfile read string
    :param read:
    :param rNext:
    :param pNext:
    :return:
    '''
    query = read[0].query
    if read[1] == True:
        for tup in read[2]:
            query = query[0:tup[0]] + tup[1] + query[tup[0] + 1:]
    string = ""
    string += str(read[0].query_name) + "\t" + str(read[0].flag) + '\t' + str(read[0].reference_name) + '\t'
    string += str(read[0].get_reference_positions()[0]) + "\t" + str(read[0].mapping_quality) + '\t' + str(
        read[0].cigarstring) + '\t' + str(rNext) + '\t' + str(pNext) + '\t' + str(read[0].template_length) + '\t' + str(
        query) + '\t' + str(read[0].qual)

    return string


def read_coverage_list(path):
    list_of_n = []
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
                continue
            list_of_n.append(int(row[1]))
            line_count += 1
    return list_of_n


def read_percentage_list(path):
    list_of_p = []
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
                continue
            list_of_p.append(float(row[1]))
            line_count += 1
    return list_of_p


def read_error_rate_list(path):
    list_of_p = []
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
                continue
            list_of_p.append(float(row[0]))
            line_count += 1
    return list_of_p


if __name__ == "__main__":
    initial_time = time.time()
    # normal_path = sys.argv[1]
    # tumor_path = sys.argv[2]
    # output_name = sys.argv[3]
    output_name = 'Simulator_v2.0.1'
    tumor_path = "test_sam.sam"
    normal_path = "p1_blood.sam"
    tumor_sam = pysam.AlignmentFile(tumor_path, "r")
    normal_sam = pysam.AlignmentFile(normal_path, "r")
    all_reads = {}
    list_of_p_tumor = read_percentage_list('ctPercentage_sample.csv')
    list_of_normal_error_rates = read_error_rate_list('bloodSmallPeak.csv')
    list_of_tissue_error_rates = read_error_rate_list('tissueSmallPeak.csv')
    window_size = 300
    p_tumor = 0.2
    number_of_input_normal_reads = 0
    number_of_input_tissue_reads = 0
    number_of_simulated_normal_reads = 0
    number_of_simulated_tissue_reads = 0
    total_positions = 1250103
    if len(list_of_normal_error_rates) < total_positions:
        total_positions = len(list_of_normal_error_rates) - 1
    header = normal_sam.text
    print("Reading normal file...")
    for normal_read in normal_sam.fetch():
        try:
            if len(normal_read.get_reference_positions()) == 0:
                # print("XXX")
                continue
            if not is_a_valid_read(normal_read.reference_name):
                # print("YYY")
                continue
            read = ['normal', normal_read.reference_name, normal_read.get_reference_positions()[0],
                    normal_read.get_reference_positions()[-1], normal_read]
            if normal_read.reference_name in all_reads:
                all_reads[normal_read.reference_name].append(read)
                number_of_input_normal_reads += 1
            else:
                all_reads.update({normal_read.reference_name: [read]})
                number_of_input_normal_reads += 1
        except AttributeError as e:
            # print(e)
            continue
    print("Reading tumor file...")
    for tumor_read in tumor_sam:
        try:
            if len(tumor_read.get_reference_positions()) == 0:
                continue
            if not is_a_valid_read(tumor_read.reference_name):
                continue

            read = ['tumor', tumor_read.reference_name, tumor_read.get_reference_positions()[0],
                    tumor_read.get_reference_positions()[-1], tumor_read]
            if tumor_read.reference_name in all_reads:
                all_reads[tumor_read.reference_name].append(read)
                number_of_input_tissue_reads += 1
            else:
                all_reads.update({tumor_read.reference_name: [read]})
                number_of_input_tissue_reads += 1
        except AttributeError as e:
            # print(e)
            continue
    print("number of total input reads:\t", number_of_input_tissue_reads + number_of_input_normal_reads)
    print("number of normal input reads:\t", number_of_input_normal_reads)
    print("number of tissue input reads:\t", number_of_input_tissue_reads)
    sam_file_reads = []
    print("Creating samfile...")
    for chr in all_reads:
        chr_reads = sorted(all_reads[chr], key=lambda x: x[2])
        print("\t\tSimulating chromosome:\t", chr)
        begin_pos = chr_reads[0][2]
        end_pos = chr_reads[-1][3]
        positions = sorted(random.sample(range(begin_pos, end_pos), int(total_positions / 24)))
        list_of_n = read_coverage_list('chr_cov/chr_cov/{0}_coverage.csv'.format(chr))
        # error_positions = [[72550, 0.5, 0.5], [110400, 0.5, 0.5], [142700, 0.5, 0.5], [142800, 0.5, 0.5],
        #                    [142900, 0.5, 0.5], [143100, 0.5, 0.5], [143300, 0.5, 0.5], [144100, 0.5, 0.5],
        #                    [144300, 0.5, 0.5]]
        error_positions = []
        for i in range(len(positions)):
            error_positions.append([positions[i], list_of_normal_error_rates[i], list_of_tissue_error_rates[i]])
        l = len(error_positions)
        list_of_normal_error_rates = list_of_normal_error_rates[l:]
        list_of_tissue_error_rates = list_of_tissue_error_rates[l:]
        read_in_range, normal_sim, tissue_sim = get_read_in_range(chr_reads, window_size, list_of_n, list_of_p_tumor,
                                                                  error_positions)
        number_of_simulated_normal_reads += normal_sim
        number_of_simulated_tissue_reads += tissue_sim
        for r in read_in_range:
            sam_file_reads.append(r)
    print("number of total simulated reads:\t", len(sam_file_reads))
    print("number of normal simulated reads:\t", number_of_simulated_normal_reads)
    print("number of tissue simulated reads:\t", number_of_simulated_tissue_reads)
    print("Writing...")
    with open('{0}_simulated.sam'.format(output_name), 'w') as writer:
        writer.write(header)
        for i in range(len(sam_file_reads)):
            read = sam_file_reads[i]
            if i < len(sam_file_reads) - 1:
                pNext = str(sam_file_reads[i + 1][0].get_reference_positions()[0])
            else:
                pNext = '0'
            writer.write(read_to_string(read, '=', pNext) + '\n')
    writer.close()
    # with open('{0}_simulated.fq'.format(output_name), 'w') as writer2:
    #     for read in sam_file_reads:
    #         writer2.write("@{0}\n".format(str(read[0].query_name)))
    #         writer2.write("{0}\n".format(str(read[0].query)))
    #         writer2.write("+\n")
    #         writer2.write("{0}\n".format(str(read[0].qual)))
    # writer2.close()

    print("Run time duration:\t{0}".format(time.time() - initial_time))
