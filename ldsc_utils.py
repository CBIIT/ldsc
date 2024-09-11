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
        command = f"cd {fileDir} && python /app/ldsc.py --bfile {file_chromo} --l2 --ld-wind-cm 1 --out {file_chromo}"
        result = subprocess.run(
            ['bash', '-c', command],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"An error occurred: {e.stderr}"

