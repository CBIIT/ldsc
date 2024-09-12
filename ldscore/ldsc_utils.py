import os
import glob
import subprocess



def run_ldsc_command(pop, genome_build, filename):
    fileDir = f"/data/tmp/uploads"
    print(filename)
    if filename:
        file_parts = filename.split('.')
        file_chromo = None
        for part in file_parts:
            if part.isdigit() and 1 <= int(part) <= 22:
                file_chromo = part
                break
    
    if file_chromo:
        # Find the file in the directory
        pattern = os.path.join(fileDir, f"{filename}.*")
        for file_path in glob.glob(pattern):
            extension = file_path.split('.')[-1]
            new_filename = f"{file_chromo}.{extension}"
            new_file_path = os.path.join(fileDir, new_filename)
            os.rename(file_path, new_file_path)
        
    try:
        # Run the command
        # 'cd 1kg_eur && python ../ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22'
        script_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(script_dir)
        ldsc_script_path = os.path.join(parent_dir, 'ldsc.py')
        #print(ldsc_script_path)
        command = f"cd {fileDir} && python3 {ldsc_script_path} --bfile {file_chromo} --l2 --ld-wind-cm 1 --out {file_chromo}"
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

