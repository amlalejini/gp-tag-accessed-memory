'''
Script: agg_min_correct_networks.py

For each run, grab smallest correct solution network. If run has none, report none.

'''

import argparse, os, copy, errno, csv

problem_whitelist = ["grade", "number-io", "for-loop-index", "median", "smallest", "small-or-large", "compare-string-lengths", "sum-of-squares", "string-lengths-backwards"]

def mkdir_p(path):
    """
    This is functionally equivalent to the mkdir -p [fname] bash command
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")
    parser.add_argument("dump_directory", type=str, help="Where to dump this?")
    parser.add_argument("-u", "--update", type=int, help="max update to look for solutions")

    args = parser.parse_args()

    data_directory = args.data_directory
    dump = args.dump_directory

    print("Pulling smallest network solutions from all runs in {}".format(data_directory))

    mkdir_p(dump)

    # Get a list of all runs
    runs = [d for d in os.listdir(data_directory) if d.strip("PROBLEM_").split("__")[0] in problem_whitelist]
    runs.sort()

    if args.update != None:
        update = args.update
        print("Looking for best solutions from update {} or earlier.".format(update))

        solutions_content = "treatment,run_id,problem,arg_type,num_arg_mut_rate,arg_tag_bf_mut_rate,arg_tag_rand_rate,mem_searching,register_tags_init,register_tags_evolve,register_capacity_evolve,register_duplication_rate,register_deletion_rate,solution_found,solution_length,update_found,update_first_solution_found,program,program_register_cnt,program_register_tags\n"

        for run in runs:
            print("Run: {}".format(run))
            run_dir = os.path.join(data_directory, run)
            run_id = run.split("__")[-1]
            run = "__".join(run.split("__")[:-1])
            treatment = run
            run_sols = os.path.join(run_dir, "output", "solutions.csv")

            # gimme dat run log
            run_log_fpath = os.path.join(run_dir, "run.log")
            with open(run_log_fpath, "r") as logfp:
                log_content = logfp.read()

            run_settings = {}
            settings_content = log_content.split("==============================")[2]
            for line in settings_content.split("\n"):
                if line[:3] == "set":
                    line = line.split(" ")
                    run_settings[line[1]] = line[2]
            # Extract relevant run settings!
            # problem,
            problem = run_settings["PROBLEM"]

            # arg_type,
            arg_type = "UNKNOWN"
            if run_settings["PROGRAM_ARGUMENT_MODE"] == "0":
                arg_type = "TAG"
            elif run_settings["PROGRAM_ARGUMENT_MODE"] == "1":
                arg_type = "NUMERIC"
            elif run_settings["PROGRAM_ARGUMENT_MODE"] == "2":
                arg_type = "BOTH"

            # num_arg_mut_rate,
            num_arg_mut_rate = run_settings["PROG_MUT__PER_NUMERIC_ARG_SUB"]

            # arg_tag_bf_mut_rate,
            arg_tag_bf_mut_rate = run_settings["PROG_MUT__PER_BIT_FLIP"]

            # arg_tag_rand_rate,
            arg_tag_rand_rate = run_settings["PROG_MUT__PER_TAG_RANDOMIZE"]

            # mem_searching,
            mem_searching = run_settings["PROGRAM_ARGUMENTS_TYPE_SEARCH"]

            # register_tags_init
            register_tags_init = "UNKNOWN"
            if run_settings["MEM_TAG_INIT_MODE"] == "0":
                register_tags_init = "hadamard"
            elif run_settings["MEM_TAG_INIT_MODE"] == "1":
                register_tags_init = "random"

            # register_tags_evolve,
            register_tags_evolve = run_settings["MEM_TAG_EVOLVE"]

            # register_capacity_evolve,
            register_capacity_evolve = run_settings["MEM_CAPACITY_EVOLVE"]

            # register_duplication_rate,
            register_duplication_rate = run_settings["MEM_TAG_MUT__PER_TAG_DUP"]

            # register_deletion_rate,
            register_deletion_rate = run_settings["MEM_TAG_MUT__PER_TAG_DEL"]

            # Extract solution information

            file_content = None
            with open(run_sols, "r") as fp:
                file_content = fp.read().strip().split("\n")
            header = file_content[0].split(",")
            header_lu = {header[i].strip():i for i in range(0, len(header))}
            file_content = file_content[1:]
            solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
            # Add smallest solution to smallest solution doc
            solution_found = False
            update_first_solution_found = "NONE"
            solution_length = "NONE"
            program = "NONE"
            program_register_cnt = "NONE"
            program_register_tags = "NONE"
            if len(solutions) > 0:
                solution_found = True
                update_first_solution_found = solutions[0][header_lu["update"]]
                solution_length = solutions[0][header_lu["program_len"]]
                program = solutions[0][header_lu["program"]]
                program_register_cnt = solutions[0][header_lu["program_register_cnt"]]
                program_register_tags = solutions[0][header_lu["program_register_tags"]]

            #"treatment,run_id,problem,arg_type,num_arg_mut_rate,arg_tag_bf_mut_rate,arg_tag_rand_rate,mem_searching,register_tags_init,register_tags_evolve,register_capacity_evolve,register_duplication_rate,register_deletion_rate,solution_found,solution_length,update_found,update_first_solution_found,program,program_register_cnt,program_register_tags\n"
            solutions_content += ",".join(map(str,[treatment,run_id,problem,arg_type,num_arg_mut_rate,arg_tag_bf_mut_rate,arg_tag_rand_rate,mem_searching,register_tags_init,register_tags_evolve,register_capacity_evolve,register_duplication_rate,register_deletion_rate,solution_found,solution_length,solution_found,update_first_solution_found,'"{}"'.format(program),program_register_cnt,'"{}"'.format(program_register_tags)])) + "\n"
        with open(os.path.join(dump, "solutions__update_{}.csv".format(update)), "w") as fp:
            fp.write(solutions_content)

if __name__ == "__main__":
    main()