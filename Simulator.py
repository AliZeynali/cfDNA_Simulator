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


def get_read_in_range(all_reads, window_size, list_of_error_positions):
    '''

    :param all_reads: list of [ type (tumor or normal ), chromosome, begin_pos, end_pos, raw_read_string]
    :param begin: first position in this chromosome
    :param window_size:
    :param list_of_error_positions: list of pos of errors
    :return:
    '''
    p_shape, p_scale = 0.7532623, 3.4025548  # For the p_tumor distribution
    n_size, n_prob = 138.7844083, 0.2221451
    cfDNA_mu, cfDNA_sigma = -1.650449, 0.3348813
    tissue_mu, tissue_sigma = -2.161420, 0.4807095
    blood_mu, blood_sigma = -1.903177, 0.3766689
    reads = []
    normal_simulated = 0
    tissue_simulated = 0
    new_position_error = []
    while len(all_reads) > 0:
        p_tumor_val = -1
        while (p_tumor_val <= 0 or p_tumor_val >= 1):
            p_tumor_val = numpy.random.gamma(p_shape, p_scale, 1)[0]
        p_normal = 1 - p_tumor_val
        i = int(all_reads[0][2] / window_size)
        n = 0
        while n < 6:
            n = numpy.random.negative_binomial(n_size, n_prob, 1)[0]
        begin_pos = i * window_size
        end_pos = begin_pos + window_size
        window_read_tumor = []
        window_read_normal = []
        if len(list_of_error_positions) > 0 and begin_pos <= list_of_error_positions[0] <= end_pos:
            list_val = list_of_error_positions.pop(0)
            p_blood_error = -1
            p_tissue_error = -1
            while (p_blood_error <= 0 or p_blood_error >= 1):
                cfdna_error = 10 ** (numpy.random.normal(cfDNA_mu, cfDNA_sigma, 1)[0])
                blood_error = 10 ** (numpy.random.normal(blood_mu, blood_sigma, 1)[0])
                p_blood_error = cfdna_error - blood_error
            while (p_tissue_error <= 0 or p_tissue_error >= 1):
                cfdna_error = 10 ** (numpy.random.normal(cfDNA_mu, cfDNA_sigma, 1)[0])
                tissue_error = 10 ** (numpy.random.normal(tissue_mu, tissue_sigma, 1)[0])
                p_tissue_error = cfdna_error - tissue_error
            new_position_error.append([list_val, p_normal * p_blood_error + p_tumor_val * p_tissue_error])
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
    query = read[0].seq
    if read[1] == True:
        for tup in read[2]:
            query = query[0:tup[0]] + tup[1] + query[tup[0] + 1:]
    qual = str(read[0].qual)
    if len(qual) == len(query) + 3:
        qual = qual[2:-1]
    elif len(qual) != len(query):
        qual = ""
        for _ in range(len(query)):
            qual += "J"

    string = ""
    string += str(read[0].query_name) + "\t" + str(read[0].flag) + '\t' + str(read[0].reference_name) + '\t'
    string += str(read[0].get_reference_positions()[0]) + "\t" + str(read[0].mapping_quality) + '\t' + str(
        read[0].cigarstring) + '\t' + str(rNext) + '\t' + str(pNext) + '\t' + str(read[0].template_length) + '\t' + str(
        query) + '\t' + qual
    # if len(query) != len(str(read[0].qual)[2:-1]):
    #     print("Hi")

    return string


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
    version = "v2.2.1"
    initial_time = time.time()
    normal_path = sys.argv[1]
    tumor_path = sys.argv[2]
    # output_name = sys.argv[3]
    output_name = 'Simulator' + "_" + version
    # tumor_path = "test_sam.sam"
    # normal_path = "p1_blood.sam"
    print("Running Simulator version:\t" + version)
    tumor_sam = pysam.AlignmentFile(tumor_path, "r")
    normal_sam = pysam.AlignmentFile(normal_path, "r")
    all_reads = {}
    window_size = 300
    p_tumor = 0.2
    number_of_input_normal_reads = 0
    number_of_input_tissue_reads = 0
    number_of_simulated_normal_reads = 0
    number_of_simulated_tissue_reads = 0
    total_positions = 1250103
    header = normal_sam.text
    error_poses = []
    print("Reading normal file...")
    for normal_read in normal_sam.fetch():
        try:
            if len(normal_read.get_reference_positions()) == 0:
                continue
            if not is_a_valid_read(normal_read.reference_name):
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
    del (normal_sam)

    sam_file_reads = []
    print("Simulating by normal reads...")
    for chr in all_reads:
        chr_reads = sorted(all_reads[chr], key=lambda x: x[2])
        print("\t\tSimulating chromosome:\t", chr, " in normal reads")
        begin_pos = chr_reads[0][2]
        end_pos = chr_reads[-1][3]
        positions = sorted(random.sample(range(begin_pos, end_pos), int(total_positions / 48)))
        for p in positions:
            error_poses.append([chr,p])
        error_positions = []
        for i in range(len(positions)):
            error_positions.append(positions[i])
        read_in_range, normal_sim, tissue_sim = get_read_in_range(chr_reads, window_size,
                                                                  error_positions)
        number_of_simulated_normal_reads += normal_sim
        number_of_simulated_tissue_reads += tissue_sim
        for r in read_in_range:
            sam_file_reads.append(r)

    del (all_reads)
    all_reads = {}
    print("Reading tumor file...")
    for tumor_read in tumor_sam:
        try:
            if len(tumor_read.get_reference_positions()) == 0:
                continue
            if not is_a_valid_read(tumor_read.reference_name):
                continue

            read = ['tumor', tumor_read.reference_name, tumor_read.get_reference_positions()[0],
                    tumor_read.get_reference_positions()[-1], tumor_read]
            if tumor_read.cigarstring != '150M':
                pass
            if tumor_read.reference_name in all_reads:
                all_reads[tumor_read.reference_name].append(read)
                number_of_input_tissue_reads += 1
            else:
                all_reads.update({tumor_read.reference_name: [read]})
                number_of_input_tissue_reads += 1
        except AttributeError as e:
            # print(e)
            continue
    del (tumor_sam)
    print("Simulating by tissue reads...")
    for chr in all_reads:
        chr_reads = sorted(all_reads[chr], key=lambda x: x[2])
        print("\t\tSimulating chromosome:\t", chr, " in tumor reads")
        begin_pos = chr_reads[0][2]
        end_pos = chr_reads[-1][3]
        positions = sorted(random.sample(range(begin_pos, end_pos), int(total_positions / 48)))
        for p in positions:
            error_poses.append([chr,p])
        error_positions = []
        for i in range(len(positions)):
            error_positions.append(positions[i])
        read_in_range, normal_sim, tissue_sim = get_read_in_range(chr_reads, window_size,
                                                                  error_positions)
        number_of_simulated_normal_reads += normal_sim
        number_of_simulated_tissue_reads += tissue_sim
        for r in read_in_range:
            sam_file_reads.append(r)

    print("\nnumber of total input reads:\t", number_of_input_tissue_reads + number_of_input_normal_reads)
    print("number of normal input reads:\t", number_of_input_normal_reads)
    print("number of tissue input reads:\t", number_of_input_tissue_reads)

    print("number of total simulated reads:\t", len(sam_file_reads))
    print("number of normal simulated reads:\t", number_of_simulated_normal_reads)
    print("number of tissue simulated reads:\t", number_of_simulated_tissue_reads)
    print("Writing to file...")
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
    with open('Error_positions.txt', 'w') as writer2:
        for val in error_poses:
            writer2.write(str(val[0]) + "\t" + str(val[1]) + "\n")

    print("Run time duration:\t{0}".format(time.time() - initial_time))
