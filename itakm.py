import re
import sys


"""
 Function: parse hmmscan result 

 Input: hmmscan result file name

 Output1: detail hits info, format is below:
	GeneID      PfamID		GA	    Evalue
	AT2G34620.1	PF02536.7	242.2	4.5e-72

 Output2: alignment of hits into, format is below:
	1. GeneID      -- AT1G01140.1
	2. PfamID      -- PF00069.18
	3. Query Start -- 19
	4. Query End   -- 274
	5. Hit Start   -- 1
	6. Hit End     -- 260
	7. Query Seq   -- YEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVF....
	8. Alignment   -- ye++++lG+Gsf+kV  ak+  tg++ A+Kil++e+  ....
	9. Hit Seq     -- yelleklGsGsfGkVykakkkktgkkvAvKilkkeeek....
	10. GA Score   -- 241.4
	11. Evalue     -- 6.9e-72 
	12. Description-- Protein kinase domain
	13. Query Len  -- 448

 Output3: for some hits, the best hits domain may not meet the requirement of significant
          but multi-hits for one domain actually meet the GA cutoff. That why the GA score
          has two values. One for overall, another for best domain
	GeneID      PfamID      GA      Evalue
    AT2G34620.1 PF02536.7   242.2   4.5e-72
=cut

"""
def parse_hmmscan_result(hmm_result):
    result_out1 = ""
    result_out2 = ""
    result_out3 = ""
    one_result = ""
    align_detail = ""
    hits = ""

    try:
        with open(hmm_result, 'r') as rfh:
            for line in rfh:
                if not line.startswith("#"):
                    if line.startswith("//"):
                        hits_content = one_result.split(">>")

                        hit_head = hits_content[0].split("\n")
                        query_name = None
                        query_length = None
                        jumper = 0

                        """
                        # parse domain info
                        1. short info for every hsp	
                        $hits = query_id \t domain_id \t evalue \t score desc \n #
                        
                        # parse hit head content like below	
                        # Match a line like this
                        E-value  score  bias    E-value  score  bias    exp  N   Sequence Description
                        -------  -----  -----   -------  ------ -----   ---- --  -------- -----------
                        4e-83   285.8  10.0    5.3e-83  285.5  7.0     1.1  1   Q14SN3.1 Q14SN3_9HEPC Polyprotein (Fragment).
                        #######################################################################################################
                        # previous hmmer3 parse function need this part for best one domain.
                        # New version just need Query name, length and no hit information in this part.
                        #######################################################################################################
                        """
                        for hit_head_line in hit_head:
                            match_query_m = re.search(r"^Query:\s+(\S+)\s+\[M\=(\d+)\]", hit_head_line)
                            match_query_l = re.search(r"^Query:\s+(\S+)\s+\[L\=(\d+)\]", hit_head_line)
                            match_no_hits = re.search(r"No hits detected that satisfy reporting thresholds", hit_head_line)

                            if match_query_m:
                                query_name, query_length = match_query_m.groups()
                            elif match_query_l:
                                query_name, query_length = match_query_l.groups()
                            elif match_no_hits:
                                jumper = 1
                            else:
                                b = hit_head_line.split(None, 11)
                                if len(b) > 2 and re.match(r"^\d+\.\d+", b[1]):
                                    result_out3 += "{}\t{}\t{}\t{}\n".format(query_name, b[8], b[1], b[0])

                        """ parse hsp part of hits	"""
                        if jumper != 1:
                            for ih in range(1, len(hits_content)):
                                one_hit = ">>" + hits_content[ih]
                                hsp_info, hsp_detail = parse_align(one_hit, query_name, query_length)
                                result_out1 += hsp_info
                                result_out2 += hsp_detail

                        one_result = ""
                    else:
                        one_result += line
    except IOError as e:
        print(f"Cannot open hmmscan result file: {hmm_result} {e}")

    return result_out1, result_out2, result_out3

# Function to parse alignment details
def parse_align(hsp_info, query_name, query_length):
    output1 = ""
    output2 = ""
    #print(hsp_info)
    hsp_lines = hsp_info.strip().split('\n')
    info1 = {}
    hit_id = hit_desc = ""
    
    for line in hsp_lines:
        if line.startswith('>> '):
            parts = line.split(None, 2)
            hit_id = parts[1]
            hit_desc = parts[2]
        elif re.match(r'^\s+\d+\s+\W\s+', line):
            info1[int(line.split()[0])] = line
    
    domains = hsp_info.split('== domain')
    #print(hsp_info)
    #print("dddd")
    #print(domains[1])
    #sys.exit(1)
    
    for j in range(1, len(domains)):
        domain_lines = domains[j].strip().split('\n')
        info = info1[j].split()

        #print(info)
        #sys.exit(1)
        
        query_string = hit_string = match_string = ""
        hsp_length = align_pos = match_start = None
        
        for line in domain_lines:
            if re.match(rf'^\s+{re.escape(hit_id)}\s+\d+\s+\S+\s+\d+', line):
                fff = line.split()
                hit_string = fff[3]
                hsp_length = len(hit_string)
                align_pos = line.find(hit_string)
                match_start = True
                if align_pos < 0:
                    raise ValueError(f"Error! Align start position is negative: {align_pos}\n{line}\n{hit_string}\n")
            elif match_start:
                match_string = line[align_pos:align_pos + hsp_length]
                match_start = False
            elif re.match(rf'^\s+{re.escape(query_name)}', line):
                qqq = line.split()
                query_string = qqq[3]
        
        output1 += f"{query_name}\t{hit_id}\t{info[2]}\t{info[5]}\n"
        output2 += f"{query_name}\t{hit_id}\t{info[9]}\t{info[10]}\t{info[6]}\t{info[7]}\t{query_string}\t{match_string}\t{hit_string}\t{info[2]}\t{info[5]}\t{hit_desc}\t{query_length}\n"
        
    return output1, output2

"""
def parse_align(hsp_info, query_name, query_length):
    output1 = ""
    output2 = ""

    hsp_line = hsp_info.split("\n")
    hit_id = None
    hit_desc = None
    info1 = {}

    for i, line in enumerate(hsp_line):
        if line.startswith(">>"):
            aa = line.split(None, 2)
            hit_id, hit_desc = aa[1], aa[2]
        elif re.match(r"^\s+(\d+)\s+\W\s+", line):
            info1[int(line.split(None, 1)[0])] = line

    for j in range(1, len(hsp_line)):
        if j in info1:
            domain_lines = info1[j].split("\n")
            info = domain_lines[0].split()
            output1 += "{}\t{}\t{}\t{}\n".format(query_name, hit_id, info[3], info[6])
            output2 += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                query_name, hit_id, info[10], info[11], info[7], info[8], query_string, match_string, hit_string, info[3], info[6], hit_desc, query_length
            )

    return output1, output2
"""


def parse_format_result(in_file, ga_score_hash):
    out_hash = {}
    hsp_hit_id = {}
    hsp_detail = {}

    in_file_lines = in_file.strip().split('\n')
    uid = 0
    len_padding = len(str(len(in_file_lines)))

    for line in in_file_lines:
        fmm = line.strip().split('\t')
        domain_id = fmm[1].split('.')[0]  # Remove anything after the dot

        valued = False
        ga_score = ga_score_hash.get(domain_id)

        if ga_score is not None:
            if float(fmm[9]) >= ga_score:
                valued = True
        else:
            if float(fmm[10]) <= 1e-3:
                valued = True

        if valued:
            uid += 1
            zero_padding = str(uid).zfill(len_padding)

            hsp_detail[zero_padding] = line
            hsp_hit_id[zero_padding] = domain_id

            if fmm[0] in out_hash:
                out_hash[fmm[0]] += "\t" + zero_padding
            else:
                out_hash[fmm[0]] = zero_padding

    return out_hash, hsp_hit_id, hsp_detail
