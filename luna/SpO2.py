import mne
import os
import glob
import numpy as np
from concurrent.futures import ThreadPoolExecutor

directory_path = '/data/nsrr/working/pats/dist'
edf_files = glob.glob(os.path.join(directory_path, '*.edf'))
output_file_path = '/data/nsrr/working/pats/scripts/output_summary2.txt'
match_strings = ["N1", "N2", "N3", "R"]
epoch_length = 30

def process_file(edf_file):    
    print(f"Processing edf file: {edf_file}")
    subject_name = os.path.basename(edf_file).split('.')[0]

    # Read the EDF file
    channels_to_include = ['SpO2']
    raw_data = mne.io.read_raw_edf(edf_file, include=channels_to_include, preload=True)

    # Select sleep indices from eannot files
    eannot_file = os.path.join('/data/nsrr/working/pats/eannot/', f'{subject_name}.eannot')
    print(f"Processing eannot file: {eannot_file}")
    string_array = np.genfromtxt(eannot_file, dtype=str, delimiter=',')
    mask = np.isin(string_array, match_strings)
    indices = np.where(mask)[0]

    if 'SpO2' in raw_data.ch_names:

        events = mne.make_fixed_length_events(raw_data, id=1, start=0, duration=epoch_length)
        epochs = mne.Epochs(raw_data, events, event_id=1, tmin=0, tmax=epoch_length - 1 / raw_data.info['sfreq'], baseline=None, preload=True)
        selected_epochs = epochs[indices]
        data = selected_epochs.get_data(copy=False)
        data_reshaped = np.concatenate(data, axis=-1)
        info = mne.create_info(ch_names=epochs.ch_names, sfreq=epochs.info['sfreq'], ch_types=epochs.get_channel_types())
        raw_from_epochs = mne.io.RawArray(data_reshaped, info)
    
        # To access the data
        spo2_data = raw_from_epochs.get_data()[0]
        samples_below_75 = sum(spo2_data < 75)
        sample_rate = raw_from_epochs.info['sfreq']
        total_seconds_below_75 = samples_below_75 / sample_rate
        total_minutes_below_75 = total_seconds_below_75 / 60
        print(f"Total minutes with SpO2 < 75: {total_minutes_below_75:.2f} minutes")

        total_samples = raw_from_epochs.get_data().shape[1]
        total_seconds = total_samples / sample_rate
        total_minutes = total_seconds / 60
        print(f"Total recording length of SpO2 in minutes: {total_minutes:.2f}")
        print("\n")
        return f"{subject_name}\t{total_minutes:.2f}\t{total_minutes_below_75:.2f}\n"
        
    else:
        print("SpO2 channel not found in this subject.")
        return f"{subject_name}\tN/A\n"

        
# Main function to run the parallel processing
def main():
    # Find all EDF files in the directory
    edf_files = glob.glob(os.path.join(directory_path, '*.edf'))

    # Open the output file for writing
    with open(output_file_path, 'w') as output_file:
        output_file.write("ID\tTotal_Minutes\tSpO2 < 75\n") # Write the header

        # Setup the ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=8) as executor:
            # Map the executor to the processing function and the list of files
            results = executor.map(process_file, edf_files)

            # Write results to file as they become available
            for result in results:
                output_file.write(result)

    print("Processing completed.")

# Run the main function
if __name__ == '__main__':
    main()
