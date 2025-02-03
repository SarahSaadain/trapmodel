from common_trapmodel_scripts import *
from trapmodel_pipeline import * 


def run_modified_pipeline():

    print_info("Running trapmodel pipeline ...")

    processed_folder = get_folder_path_processed()
    if os.path.exists(processed_folder):
        print(f"Processed folder exists -> ABORT")
        return

    # Specify the folder paths
    processed_LopikProtrac_LopikCluster = os.path.join(get_folder_trapmodel(), 'processed_LopikProtrac_LopikCluster')
    processed_folder = get_folder_path_processed()

    # Rename the folder
    # processed_LopikProtrac_LopikCluster -> processed
    try:
        os.rename(processed_LopikProtrac_LopikCluster, processed_folder)
        print(f"Folder renamed from '{processed_LopikProtrac_LopikCluster}' to '{processed_folder}'")
    except Exception as e:
        print(f"Error: {e}")
        return
    
    #from trapmodel_pipeline
    run_pipeline()

    # reverse renaming
    try:
        os.rename(processed_folder, processed_LopikProtrac_LopikCluster)
        print(f"Folder renamed from '{processed_folder}' to '{processed_LopikProtrac_LopikCluster}'")
    except Exception as e:
        print(f"Error: {e}")
        return


def main():

    run_pipeline()

if __name__ == "__main__":
    main()