import requests
import tarfile
import subprocess
import os

def create_env_file(env_vars):
    """
    Creates a .env file with the given environment variables in the script's directory.

    :param env_vars: Dictionary of environment variables (key-value pairs)
    """
    # Directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Path to the .env file
    env_file_path = os.path.join(script_dir, '.env')

    # Convert the dictionary to .env format
    env_content = '\n'.join([f'{key}={value}' for key, value in env_vars.items()])

    # Create and write to the .env file
    with open(env_file_path, 'w') as file:
        file.write(env_content)

    print(f".env file created at {env_file_path}")

def download_and_extract(url, dest_folder):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dest_folder = os.path.join(script_dir, dest_folder)

    # Ensure the destination folder exists
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    # Download the file
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        file_name = url.split('/')[-1]
        tar_path = os.path.join(dest_folder, file_name)

        # Write the downloaded file to a new file in the destination folder
        with open(tar_path, 'wb') as file:
            file.write(response.raw.read())

        # Extract the .tar.gz file
        with tarfile.open(tar_path, 'r:gz') as tar:
            tar.extractall(path=dest_folder)

        # Remove the .tar.gz file after extraction
        os.remove(tar_path)
    else:
        print(f"Error downloading the file: HTTP {response.status_code}")

def find_pw_x_executable():
    """
    Finds the path of the 'pw.x' executable.
    
    Returns:
        The path of 'pw.x' if found, otherwise None.
    """
    try:
        # Unix/Linux/Mac command to find executable
        cmd = ["which", "pw.x"]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None


if __name__ == "__main__":
    download_and_extract('https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBEsol_efficiency.tar.gz', 'pp/efficiency')
    download_and_extract('https://archive.materialscloud.org/record/file?record_id=1732&filename=SSSP_1.3.0_PBEsol_precision.tar.gz', 'pp/precision')

    executable_path = find_pw_x_executable()
    script_dir = os.path.dirname(os.path.abspath(__file__))


    env_vars = {
        'PWSCF_COMMAND': executable_path,
        'PSEUDOPOTENTIALS': os.path.join(script_dir, 'pp/efficiency'),
    }

    create_env_file(env_vars)

