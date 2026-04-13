import os
import gzip
import shutil

def unzip_gz():
    base_dir = "FermiGBM_dataset"
    os.chdir(base_dir)
    print("Current Working Directory:", os.getcwd())

    # Loop over all subfolders
    for folder in os.listdir():
        folder_path = os.path.join(os.getcwd(), folder)

        if not os.path.isdir(folder_path):
            continue  # skip non-directories

        # Loop over files inside the folder
        for fname in os.listdir(folder_path):
            if fname.endswith(".gz"):
                gz_path = os.path.join(folder_path, fname)
                out_path = os.path.join(folder_path, fname[:-3])  # remove .gz

                print(f"Decompressing: {gz_path}")

                with gzip.open(gz_path, "rb") as f_in:
                    with open(out_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                print(f"Saved decompressed file as: {out_path}")

def remove_gz():
    base_dir = "FermiGBM_dataset"
    os.chdir(base_dir)
    print("Current Working Directory:", os.getcwd())

    # Loop over all subfolders
    for folder in os.listdir():
        folder_path = os.path.join(os.getcwd(), folder)

        if not os.path.isdir(folder_path):
            continue  # skip non-directories

        # Loop over files inside the folder
        for fname in os.listdir(folder_path):
            if fname.endswith(".gz"):
                gz_path = os.path.join(folder_path, fname)

                print(f"Removing: {gz_path}")
                os.remove(gz_path)
                print(f"Removed file: {gz_path}")

if __name__ == "__main__":
    #unzip_gz()
    remove_gz()
