'''
Script: agg_min_correct_networks.py

For each run, grab smallest correct solution network. If run has none, report none.

'''

import argparse, os, copy, errno, csv


def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")

    args = parser.parse_args()

    data_directory = args.data_directory
    
    # Get a list of all runs
    runs = [d for d in os.listdir(data_directory) if "__" in d]
    runs.sort()

    undone_content = "treatment,run_id,target_gen,final_update\n"
    cnt = 0
    total = 0
    for run in runs:
        print("Run: {}".format(run))
        run_dir = os.path.join(data_directory, run)
        run_id = run.split("__")[-1]
        run_name = "__".join(run.split("__")[:-1])
        run_log = None
        with open(os.path.join(run_dir, "run.log"), "r") as fp:
            run_log = fp.read().strip().split("\n")
        target_gen = None
        for line in run_log:
            if "set GENERATIONS" in line:
                target_gen = line.split(" ")[2]
                break
        if target_gen == None:
            print("Failed to find target generations in run log!")
            exit(-1)
        # Did this run finish?
        final_line = run_log[-1]
        finished = False
        if "Update: {};".format(target_gen) in final_line:
            finished = True
            print("  ==> Finished!")
        else:
            finished = False
            print("  ==> Not Finished!")
            cnt+=1
        final_update = final_line.split(",")[0].split(" ")[-1]
        
        if not finished:
            undone_content += ",".join([run_name, run_id, target_gen, final_update]) + "\n"
            
        total += 1

    with open("undone.csv", "w") as fp:
        fp.write(undone_content)
    print ("Runs not done: " + str(cnt))
        

if __name__ == "__main__":
    main()