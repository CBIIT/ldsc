import os
import glob
import subprocess
import shutil  



def run_ldsc_command(pop, genome_build, filename,ldwindow,windUnit,isExample,reference):
    fileDir = f"/data/tmp/uploads/"
    if isinstance(isExample, str):
            isExample = isExample.lower() == 'true'
  
    
    #if isExample:
    #    fileDir =  "/data/ldscore"
    ldwindow_value = 1  # Example value, replace with actual value
    # Check if ldwindow is an integer greater than 0, if not set it to 1
    try:
        ldwindow_value = int(ldwindow)
        if ldwindow_value <= 0:
            ldwindow_value = 1
    except ValueError:
        ldwindow_value = 1

    windFlag = '--ld-wind-cm'
    if windUnit == 'cm':
        windFlag = "--ld-wind-cm"
    elif windUnit == 'kb':
        windFlag = "--ld-wind-kb"

    if filename:
        #file_parts = filename.split('.')
        file_chromo = filename
        # for part in file_parts:
        #     if part.isdigit() and 1 <= int(part) <= 22:
        #         file_chromo = part
        #         break
    
    # if file_chromo:
    #     # Find the file in the directory
    #     pattern = os.path.join(fileDir, f"{filename}.*")
    #     for file_path in glob.glob(pattern):
    #         extension = file_path.split('.')[-1]
    #         new_filename = f"{file_chromo}.{extension}"
    #         new_file_path = os.path.join(fileDir, new_filename)
    #         #os.rename(file_path, new_file_path)
    #         shutil.copy(file_path, new_file_path)  # Copy the file instead of renaming it
    
    file_chr=file_chromo
    if isExample:
        file_chr =  "/data/ldscore/"+file_chromo 
        fileDir = f"/data/tmp/uploads/{reference}/"  
    else:
        fileDir = f"/data/tmp/uploads/{reference}/"   
    try:
        # Run the command
        # 'cd 1kg_eur && python ../ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22'
        parent_dir = '/usr/local/bin/'
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print(fileDir)
        command = f"cd {fileDir} && python3 {ldsc_script_path} --bfile {file_chr} --l2 {windFlag} {ldwindow_value}  --out {file_chromo}"
        result = subprocess.run(
            ['bash', '-c', command],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"An error occurred: {e.stderr}"
    

#def main():
    # Example parameters for testing
#    pop = 'example_pop'
#    genome_build = 'example_genome_build'
#    filename = 'example_filename.22.txt'
    
    # Call the function and print the result
#    result = run_ldsc_command(pop, genome_build, filename)
#    print(result)

#if __name__ == "__main__":
#    main()

def run_herit_command(sumstats_file, fileDir, ld_scores_dir,isExample):
    fallExampleDir = f"/data/ldscore"
    w_hm3_snplist = f"/data/ldscore/w_hm3.snplist"
    if isinstance(isExample, str):
        isExample = isExample.lower() == 'true'
    errormsg = ""
    try:
        parent_dir = '/usr/local/bin/'
        munge_sumstat_script_path = os.path.join(parent_dir, 'munge_sumstats.py')
        # Generate the output filename based on the input summary statistics file
        base_name = os.path.splitext(os.path.basename(sumstats_file))[0]
        out_file = f"{base_name}.sumstats.gz"
       
                # If isExample is True, use the fallbackDir
        if isExample:
            sumstats_path = os.path.join(fallExampleDir, sumstats_file)
        else:
            sumstats_path = sumstats_file
        print("First command ################:",sumstats_path, isExample)
                # Ensure ld_scores_dir is in lowercase
        ld_scores_dir = ld_scores_dir.lower()
        ld_scores_dir = fallExampleDir+"/"+ld_scores_dir
        # Ensure ld_scores_dir has a trailing slash
        if not ld_scores_dir.endswith('/'):
            ld_scores_dir += '/'
        # First command
        command1 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path} --merge-alleles {w_hm3_snplist}  --out {base_name}"
        
        #command1 = f"python ../munge_sumstats.py --sumstats {sumstats_path} --merge-alleles ../testData/w_hm3.snplist  --out {base_name}"
      
        result1 = subprocess.run( ['bash', '-c', command1], check=True, capture_output=True, text=True)
       
        # command1 = [
        #     'python', '../munge_sumstats.py',
        #     '--sumstats', sumstats_file,
        #     '--merge-alleles', '../testData/w_hm3.snplist',
        #     '--a1', 'ALT',
        #     '--a2', 'REF',
        #     '--chunksize', '500000',
        #     '--out', base_name
        # ]
        #result1 = subprocess.run( command1, check=True, capture_output=True, text=True)
       
        print("First command output:", result1.stdout)
        #print("First command error (if any):", result1.stderr)

        # Second command
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print("Second command ################:",ld_scores_dir)
        command2 = f"cd {fileDir} && python3 {ldsc_script_path} --h2 {out_file} --ref-ld-chr {ld_scores_dir} --w-ld-chr {ld_scores_dir} --out {base_name}"
       
        #result2 = subprocess.run( ['bash', '-c', command2], check=True, capture_output=True, text=True)
        try:
            result2 = subprocess.run(['bash', '-c', command2], check=True, capture_output=True, text=True)
            print("Second command output:", result2.stdout, result2.stderr)
            separator = "\n--------\n"
            return result1.stdout + separator + result2.stdout
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the second command: {e}")
            return f"{result1.stdout} An error occurred while running the second command: {e.stderr}"

        # command2 = [
        #     'python', '../ldsc.py',
        #     '--h2', "test2.sumstats.gz",
        #     '--ref-ld-chr', ld_scores_dir,
        #     '--w-ld-chr', ld_scores_dir,
        #     '--out', base_name
        # ]
        # result2 = subprocess.run( command2, check=True, capture_output=True, text=True)
       
        #print("Second command output:", result2.stdout)
        #separator = "\n--------\n"
        #return result1.stdout + separator + errormsg
        #print("Second command error (if any):", result2.stderr)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
        print(f"Command output: {e.output}")
        print(f"Command stderr: {e.stderr}")
        return f"An error occurred while running the command: {e.stderr}"

def run_correlation_command(sumstats_file, sumstats_file2, fileDir, ld_scores_dir,isExample):
    fallExampleDir = f"/data/ldscore"
    w_hm3_snplist = f"/data/ldscore/w_hm3.snplist"
    if isinstance(isExample, str):
        isExample = isExample.lower() == 'true'
    errormsg = ""
    try:
        parent_dir = '/usr/local/bin/'
        munge_sumstat_script_path = os.path.join(parent_dir, 'munge_sumstats.py')
        # Generate the output filename based on the input summary statistics file
        base_name = os.path.splitext(os.path.basename(sumstats_file))[0]
        out_file = f"{base_name}.sumstats.gz"
        base_name2 = os.path.splitext(os.path.basename(sumstats_file2))[0]
        out_file2 = f"{base_name2}.sumstats.gz"
       
                # If isExample is True, use the fallbackDir
        if isExample:
            sumstats_path = os.path.join(fallExampleDir, sumstats_file)
            sumstats_path2 = os.path.join(fallExampleDir, sumstats_file2)
        else:
            sumstats_path = sumstats_file
            sumstats_path2 = sumstats_file2
        print("First command ################:",sumstats_path, isExample)
                # Ensure ld_scores_dir is in lowercase
        ld_scores_dir = ld_scores_dir.lower()
        ld_scores_dir = fallExampleDir+"/"+ld_scores_dir
        # Ensure ld_scores_dir has a trailing slash
        if not ld_scores_dir.endswith('/'):
            ld_scores_dir += '/'
        # First command
        command1 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path} --merge-alleles {w_hm3_snplist}  --out {base_name}"
        command12 = f"cd {fileDir} && python3 {munge_sumstat_script_path} --sumstats {sumstats_path2} --merge-alleles {w_hm3_snplist}  --out {base_name2}"
        
        #command1 = f"python ../munge_sumstats.py --sumstats {sumstats_path} --merge-alleles ../testData/w_hm3.snplist  --out {base_name}"
      
        result1 = subprocess.run( ['bash', '-c', command1], check=True, capture_output=True, text=True)
        result12 = subprocess.run( ['bash', '-c', command12], check=True, capture_output=True, text=True)
      
        print("First command output:", result1.stdout)
        print("First command output:", result12.stdout)
        #print("First command error (if any):", result1.stderr)

        # Second command
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        print("Second command ################:",ld_scores_dir)
        command2 = f"cd {fileDir} && python3 {ldsc_script_path} --rg {out_file},{out_file2} --ref-ld-chr {ld_scores_dir} --w-ld-chr {ld_scores_dir} --out {base_name}"
       
        #result2 = subprocess.run( ['bash', '-c', command2], check=True, capture_output=True, text=True)

        try:
            result2 = subprocess.run(['bash', '-c', command2], check=True, capture_output=True, text=True)
            print("Second command output:", result2.stdout)
            separator = "\n--------\n"
            return result1.stdout + separator + result12.stdout + separator + result2.stdout
            #return result2.stdout
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the second command: {e}")        
            print(f"Command stderr: {e}")
            return f"An error occurred while running the second command: {e}"

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
        return f"An error occurred while running the command: {e}"


# Example usage
if __name__ == "__main__":
    user_input_sumstats = '../testData/sample/BBJ_HDLC.txt'  # Replace with actual user input
    user_input_ld_scores = '../testData/seu/'  # Replace with actual user input
    #run_herit_command(user_input_sumstats, user_input_ld_scores,False)
    #run_herit_command(user_input_sumstats, user_input_ld_scores,False)