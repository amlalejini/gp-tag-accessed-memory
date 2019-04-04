'''
Script: agg_min_correct_networks.py

For each run, grab smallest correct solution network. If run has none, report none.

'''

import argparse, os, copy, errno, csv

problem_whitelist = ["grade", "number-io", "for-loop-index", "median", "smallest", "small-or-large", "compare-string-lengths", "sum-of-squares"]

arg_types = {
    "ARGS_BOTH": "Both",
    "ARGS_NUM_ONLY": "Numeric",
    "ARGS_NUM": "Numeric",
    "ARGS_TAG_ONLY": "Tag-based",
    "ARGS_TAG_BF": "Tag-BitFlips"
}

mut_rates = {
    "MUT_005": "0.005",
    "MUT_5": "0.5",
    "MUT_1": "0.1",
    "MUT_01": "0.01",
    "MUT_001": "0.001",
    "MUT_0001": "0.0001",
    "MUT_00001": "0.00001"
}

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

        solutions_content = "treatment,run_id,problem,arg_type,arg_mut_rate,mem_searching,solution_found,solution_length,update_found,update_first_solution_found,program\n"

        for run in runs:
            print("Run: {}".format(run))
            run_dir = os.path.join(data_directory, run)
            run_id = run.split("__")[-1]
            run = "__".join(run.split("__")[:-1])
            treatment = run
            run_sols = os.path.join(run_dir, "output", "solutions.csv")

            problem = run.strip("PROBLEM_").split("__")[0]

            arg_type = None
            for thing in arg_types:
                if thing in treatment: arg_type = arg_types[thing]
            if arg_type == None:
                print("Unrecognized arg type! Exiting.")
                exit()

            arg_mut_rate = None
            for thing in mut_rates:
                if thing in treatment: arg_mut_rate = mut_rates[thing]
            if arg_mut_rate == None:
                print("Unrecognized arg mut rate! Exiting.")
                exit()

            mem_searching = "0" if "MEM_SEARCH_0" in treatment else "1"

            file_content = None
            with open(run_sols, "r") as fp:
                file_content = fp.read().strip().split("\n")

            header = file_content[0].split(",")
            header_lu = {header[i].strip():i for i in range(0, len(header))}
            file_content = file_content[1:]

            solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
            # Add smallest solution to smallest solution doc
            min_program = None
            sol_found = False
            if len(solutions) > 0:
                # Find the smallest program
                for i in range(0, len(solutions)):
                    sol_update = int(solutions[i][header_lu["update"]])
                    if sol_update > update: continue
                    if min_program == None:
                        min_program = i
                        sol_found = True
                    elif float(solutions[i][header_lu["program_len"]]) < float(solutions[min_program][header_lu["program_len"]]):
                        min_program = i
                        sol_found = True

            if sol_found:
                # Record timing info about first solution
                update_first_sol = solutions[0][header_lu["update"]]
                # Record info about smallest solution
                min_sol = solutions[min_program]
                program_len = min_sol[header_lu["program_len"]]
                update_found = min_sol[header_lu["update"]]
                program = min_sol[header_lu["program"]]
            else:
                update_first_sol = "NONE"
                program_len = "NONE"
                update_found = "NONE"
                program = "NONE"
            # "treatment,run_id,problem,uses_cohorts,solution_found,solution_length,update_found,program\n"
            solutions_content += ",".join(map(str,[treatment, run_id, problem, arg_type, arg_mut_rate, mem_searching, sol_found, program_len, update_found, update_first_sol, '"{}"'.format(program)])) + "\n"
        with open(os.path.join(dump, "min_programs__update_{}.csv".format(update)), "w") as fp:
            fp.write(solutions_content)

if __name__ == "__main__":
    main()