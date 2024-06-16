import os
import tkinter as tk
from tkinter import filedialog
import subprocess
import time
import threading
import paramiko
import sys
import glob
import datetime
from ROOT import TFile, TCanvas, TH2F, TDatime, gApplication, gStyle, gPad
import numpy as np
import re
import concurrent.futures
import time
from iqtools import *
from datetime import datetime
import toml
from PyQt5.QtWidgets import (QApplication, QWidget, QLabel, QLineEdit, QPushButton,
                             QVBoxLayout, QHBoxLayout, QMessageBox, QComboBox, QFormLayout)
from PyQt5.QtCore import QThread, pyqtSignal
def filter_files_by_time(directory, prefix, extension, start_time_str, end_time_str, time_format="%Y_%m_%d_%H_%M_%S"):
    """
    Filters ROOT files in a directory by a time range specified in their filenames.
    
    :param directory: Path to the directory containing ROOT files.
    :param start_time_str: Start time as a string.
    :param end_time_str: End time as a string.
    :param time_format: The datetime format used in the filenames.
    :return: A list of file paths that fall within the specified time range.
    """
    #print("chenrj filter_files_by_time ..."," start_time_str ",start_time_str," end_time_str ",end_time_str)
    start_time = datetime.datetime.strptime(start_time_str, time_format)
    end_time   = datetime.datetime.strptime(end_time_str,   time_format)
    #print("chenrj filter_files_by_time ..."," start_time ",start_time," end_time ",end_time)
    # Find all ROOT files in the specified directory
    #root_files = glob.glob(os.path.join(directory, extension))
    root_files = sorted(glob.glob(os.path.join(directory, extension)))
    #print("directory ", directory," extension ",extension," root_files ",root_files)
    #for file in root_files:
    #    print("root_file ",file)
    filtered_files = []
    # Use a list to store (file path, datetime) tuples
    files_with_times = []
    for file in root_files:
        # Simplify the extraction of the timestamp from the filename
        cleaned_extension = extension.replace("*","")
        filename_without_extension = os.path.basename(file).replace(cleaned_extension, "")

        #print("filename_without_extension ",filename_without_extension)
        
        timestamp_str = filename_without_extension.replace(prefix, "")
        
        #print("timestamp_str ",timestamp_str," time_format ",time_format, " prefix=",prefix,".")
        file_time = datetime.datetime.strptime(timestamp_str, time_format)
        
        # Check if the file's timestamp falls within the start and end times
        if start_time <= file_time <= end_time:
            #filtered_files.append(file)
            files_with_times.append((file, file_time))
            
        # Sort the list by the datetime object
        files_with_times.sort(key=lambda x: x[1])
        # Extract the sorted file paths
        filtered_files = [file for file, _ in files_with_times]
    #print("start_time ",start_time)
    #print("stop_time ",start_time)
    #print("filtered_files ",filtered_files)
    #for filename in filtered_files:
    #    print(filename)
    return filtered_files

        
def read_fft_injection_from_files(filtered_files, extension,dtime):
    """
    Reads and processes the FFT_injection objects from the filtered ROOT files.
    
    :param filtered_files: List of ROOT files to read from.
    """
    # Initialize an accumulator histogram named FFT_accumulated. Adjust the parameters
    # according to the expected binning and range of your FFT_injection histograms.
    FFT_accumulated = None
    # Check the extension of the files to determine processing method
    #print("extension = ", extension)

    # Initialize the arrays to accumulate the values
    xx, yy, zz = None, None, None
    
    if extension == "*.tiq.npz":
        for filename in filtered_files:
            print("adding :", filename)
            data = np.load(filename)
            xx0, yy0, zz0 = data['arr_0'], data['arr_1'], data['arr_2']
            # Accumulate the values
            if xx is None:
                xx, yy, zz = xx0, yy0, zz0
            else:
                zz += zz0
        FFT_injection = get_root_th2d(xx, yy, zz, name=filename, title=filename)
        if FFT_accumulated is None:
            # Clone the first histogram to create FFT_accumulated
            FFT_accumulated = FFT_injection.Clone("FFT_accumulated")
        else:
            FFT_accumulated.Add(FFT_injection)
    elif extension == "*.tiq":
        for filename in filtered_files:
            print("adding :", filename)
            iq = get_iq_object(filename)
            fs = iq.fs
            lframes = int(dtime*fs)
            
            nframes = int(iq.nsamples_total / lframes) - 1 if int(iq.nsamples_total / lframes) % 2 else int(iq.nsamples_total / lframes)
            print("fs",fs,"dtime = ",dtime,"lframes =",lframes, "nframes = ",nframes)
            iq.read(nframes=nframes, lframes=lframes)
            iq.method = 'fftw'
            xx0, yy0, zz0 = iq.get_power_spectrogram(lframes=lframes, nframes=nframes, sparse=True)
            # Accumulate the values
            if xx is None:
                xx, yy, zz = xx0, yy0, zz0
                xx = xx + iq.center
            else:
                zz += zz0
        if not filtered_files:
            raise ValueError("The filtered_files list is empty. Ensure it contains at least one file.")
            exit
        np.savez('combine_injection.npz', xx, yy, zz)
        FFT_injection_title = "from " + filtered_files[0] + "to" + filtered_files[-1]
        FFT_injection = get_root_th2d(xx, yy, zz, name="h", title=FFT_injection_title)
        if FFT_accumulated is None:
            # Clone the first histogram to create FFT_accumulated
            FFT_accumulated = FFT_injection.Clone("FFT_accumulated")
        else:
            FFT_accumulated.Add(FFT_injection)
    else:
        for file_path in filtered_files:
            f = TFile(file_path)
            FFT_injection = f.Get("FFT_injection")
            if FFT_injection:
                # Process the FFT_injection object here
                print(f"Read FFT_injection from {file_path}")
                if FFT_accumulated is None:
                    # Clone the first histogram to create FFT_accumulated
                    FFT_accumulated = FFT_injection.Clone("FFT_accumulated")
                    FFT_accumulated.SetDirectory(0)  # Remove it from the file so it doesn't get deleted with f.Close()
                else:
                    # Add the current histogram to FFT_accumulated
                    FFT_accumulated.Add(FFT_injection)
    
            f.Close()
    
    output_file = TFile("combine_injection.root", "RECREATE")
    FFT_accumulated.Write()
    output_file.Close()
        
    return FFT_accumulated
def read_injection_list(file_path):
    """
    Reads injection information from a file, returning a list of tuples.
    Each tuple contains (fileID, frameID_start, frameID_stop).
    If the file does not exist or cannot be opened, returns an empty list.
    
    Parameters:
    - file_path: str - The path to the file containing the injection data.
    
    Returns:
    - injections: list of tuples - List of tuples, each containing (fileID, frameID_start, frameID_stop).
    """
    injections = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 4:
                    injection_name, file_path, frameID,scaler_file_path = parts
                    frameID = int(frameID)
                    injections.append((injection_name, file_path, frameID))
                else:
                    print("The column number of {file_path} is not 4. Please check it.")
    except FileNotFoundError:
        print(f"Warning: The file '{file_path}' does not exist.")
    except IOError:
        print(f"Error: Could not read the file '{file_path}'.")
    return injections
                            

def load_data_from_single_file(file_path, frameID_start, frameID_stop):
    """
    Loads frequency, time, and amplitude data from separate binary files for the given frameID range.
    
    Parameters:
    - file_path: The base path to the .root file (used to construct paths to binary files).
    - frameID_start: The starting frame ID for data loading.
    - frameID_stop: The ending frame ID for data loading.
    
    Returns:
    - Tuple containing arrays for frequency data, adjusted time data, and amplitude data.
    """
    
    # Construct the paths to the binary files
    bin_fre_path  = file_path.replace(".root", ".bin_fre")
    bin_time_path = file_path.replace(".root", ".bin_time")
    bin_amp_path  = file_path.replace(".root", ".bin_amp")
    
    # Load data from the binary files
    fre  = np.fromfile(bin_fre_path, dtype=np.float64)
    time = np.fromfile(bin_time_path, dtype=np.float32)
    amp  = np.fromfile(bin_amp_path, dtype=np.float32)
    
    # Assuming that amp data needs to be reshaped according to time and frequency dimensions
    nSpecBins = len(fre)
    nTimeBins = len(time)  # Calculate the number of time bins based on amp array length and frequency bins
    amp = amp.reshape((nTimeBins, nSpecBins))
    
    # Adjust time and amp arrays based on the specified frameID range
    # This assumes that each frame corresponds to a row in the amp array and an element in the time array
    adjusted_time = time[frameID_start:frameID_stop]
    adjusted_amp = amp[frameID_start:frameID_stop, :]
    #print("chenrj adjusted_time ",adjusted_time)
    #for i in range(0,len(adjusted_time)):
    #    print(" adjusted_time i= ",i,adjusted_time[i]," adjusted_amp =",adjusted_amp[i][0])
    return fre, adjusted_time, adjusted_amp

def load_data_from_current_files_to_next_file(file_path, next_file_path, frameID_start, frameID_stop):
    """
    Loads data across files from a specified start to stop frame ID.
    """
    # Prepare paths for the binary data in the current file
    bin_fre_path = file_path.replace(".root", ".bin_fre")
    bin_time_path = file_path.replace(".root", ".bin_time")
    bin_amp_path = file_path.replace(".root", ".bin_amp")
    
    # Load data from the current file
    #print("chenrj min(frameID_stop, 4096)", min(frameID_stop, 4096))
    fre, time_current, amp_current = load_data_from_single_file(file_path, frameID_start, min(frameID_stop, 4095))
    
    #print("chenrj time_current ",time_current)
    # Initialize combined arrays with data from the current file
    combined_time = time_current
    combined_amp = amp_current
    
    # If frameID_stop exceeds the range, load data from the next file
    if frameID_stop > 4096:
        
        # Calculate the new frame indices for the next file
        next_start_frame = 0
        next_end_frame = frameID_stop - 4096
        
        # Load data from the next file
        _, time_next, amp_next = load_data_from_single_file(next_file_path, next_start_frame, next_end_frame)
        time_max = max(time_current)
        time_next = time_next + time_max
        # Combine the current and next file data
        combined_time = np.concatenate((combined_time, time_next))
        combined_amp = np.concatenate((combined_amp, amp_next))
        #print("chenrj load_data_from_current_files_to_next_file completed.")
        return fre, combined_time, combined_amp

def load_data_from_previous_files_to_current_file(previous_file_path, file_path, frameID_start, frameID_stop):
    """
    Loads data across files from a specified start to stop frame ID.
    """
    # Prepare paths for the binary data in the current file
    bin_fre_path = file_path.replace(".root", ".bin_fre")
    bin_time_path = file_path.replace(".root", ".bin_time")
    bin_amp_path = file_path.replace(".root", ".bin_amp")

    # Load data from the current file
    fre, time_current, amp_current = load_data_from_single_file(file_path, 0, frameID_stop)
    #print("chenrj frameID_stop = ", frameID_stop)
    #print("chenrj time_current = ", time_current)
    
    # If frameID_stop exceeds the range, load data from the next file
    if frameID_start <0:
        
        # Calculate the new frame indices for the next file
        previous_start_frame = frameID_start + 4096
        previous_end_frame = 4096
        #print("chenrj previous_start_frame=",previous_start_frame)
        #print("chenrj previous_end_frame = ", previous_end_frame)
        # Load data from the previous file
        _, time_previous, amp_previous = load_data_from_single_file(previous_file_path, previous_start_frame, previous_end_frame)
        
        # Initialize combined arrays with data from the current file
        time_max = max(time_previous)
        combined_time = time_previous - time_max
        combined_amp = amp_previous
        #print("chenrj combined_time ",combined_time)
        #print("chenrj time_current  ",time_current)
        # Combine the current and next file data
        combined_time = np.concatenate((combined_time, time_current))
        combined_amp = np.concatenate((combined_amp, amp_current))
        #print("chenrj combined_time ",combined_time)
        #print("chenrj load_data_from_previous_files_to_current_file completed.")
        return fre, combined_time, combined_amp
    

def adjust_and_load_frames(file_path, frameID, frameID_start, frameID_stop):
    """
    Adjusts the file path based on the frameID and specified offsets, and loads frame data.
    For frameID_stop >= 4096, data is read from the current file up to the next.
    For frameID_start <= 0, data is read from the previous file and the current file.
    
    Parameters:
    - file_path: The path to the current .root file.
    - frameID: The base frame ID from which offsets are calculated.
    - frameID_start: the start of frameID
    - frameID_stop:  the stop of frameID
    
    Returns:
    - A tuple containing frequency data, time data, and amplitude data across the needed frames.
    """
    
    # Function to generate the next or previous file path based on the current path and direction
    def get_adjacent_file_path(current_path, direction):
        # Extract the numeric part of the file name
        base_path, filename = os.path.split(current_path)
        match = re.search(r"(\d{7})\.iq\.tdms\.root$", filename)
        if not match:
            raise ValueError("File name does not match expected pattern.")
        file_num = int(match.group(1))
        # Adjust the file number based on the direction
        adjacent_file_num = file_num + direction
        # Generate the new file name and path
        new_filename = "{:07d}.iq.tdms.root".format(adjacent_file_num)
        return os.path.join(base_path, new_filename)
    
    # Initialize lists to hold combined data across files
    combined_fre, combined_time, combined_amp = [], [], []
    #print("chenrj frameID_start = ",frameID_start)
    #print("chenrj frameID_stop = ",frameID_stop)
    # Check and handle the case when frameID_stop exceeds the current file's range
    if frameID_stop >= 4096:
        next_file_path = get_adjacent_file_path(file_path, 1)  # Get the path for the next file
        # Load data from the current file (frameID_start to 4095) and the next file (0 to frameID_stop - 4096)
        # This is a placeholder for the actual data loading logic, which depends on your file format and data structure
        # Example:
        #print("chenrj file_path ",file_path)
        #print("chenrj next_file_path ",next_file_path)
        combined_fre, combined_time, combined_amp = load_data_from_current_files_to_next_file(file_path, next_file_path, frameID_start, frameID_stop)
        
        # Check and handle the case when frameID_start is less than or equal to 0
    elif frameID_start <= 0:
        previous_file_path = get_adjacent_file_path(file_path, -1)  # Get the path for the previous file
        # Load data from the previous file (4096 + frameID_start to 4095) and the current file (0 to frameID_stop)
        # This is a placeholder for the actual data loading logic
        # Example:
        #print("chenrj previous_file_path = ",previous_file_path)
        combined_fre, combined_time, combined_amp = load_data_from_previous_files_to_current_file(previous_file_path, file_path, frameID_start, frameID_stop)
        #print("chenrj combined_time = ",combined_time)
        
    else:
        # Load data normally from the current file if no adjustments are needed
        # This is a placeholder for your data loading logic when reading from a single file
        # Example:
        combined_fre, combined_time, combined_amp = load_data_from_single_file(file_path, frameID_start, frameID_stop)
        #print("chenrj combined_time = ",combined_time)
        
    return combined_fre, combined_time, combined_amp
    
def check_bin_data_files_exist(file_path):
    """
    Checks if the binary data files for frequency, time, and amplitude exist.
    
    Parameters:
    - file_path: The base path to the original .root file.
    
    Returns:
- True if all binary data files exist, False otherwise.
"""
    
    # Prepare paths for the binary data in the current file
    bin_fre_path = file_path.replace(".root", ".bin_fre")
    bin_time_path = file_path.replace(".root", ".bin_time")
    bin_amp_path = file_path.replace(".root", ".bin_amp")
    
    # Check if all files exist
    fre_exists = os.path.exists(bin_fre_path)
    time_exists = os.path.exists(bin_time_path)
    amp_exists = os.path.exists(bin_amp_path)
    # Initialize a list to keep track of missing files
    missing_files = []
        
    # Check each file and add to missing list if not found
    if not os.path.exists(bin_fre_path):
        missing_files.append(bin_fre_path)
    if not os.path.exists(bin_time_path):
        missing_files.append(bin_time_path)
    if not os.path.exists(bin_amp_path):
        missing_files.append(bin_amp_path)
        
    # Report missing files and return False if any are missing
    if missing_files:
        print("The following required data files are missing:")
        for missing_file in missing_files:
            print(missing_file)
        return False
    else:
        return True
    
def check_root_data_files_exist(injection_list_path, injection_name):
    """
    Checks if the "_full.root" data files for full spectrum exist.

    Parameters:
    - injection_name: The base path to the  "_full.root" file.

    Returns:
- True if the "_full.root" files exist, False otherwise.
"""
    # Prepare paths for the binary data in the current file
    injection_full_list_path = injection_list_path.replace("_list.txt","_full_list.txt")
    full_root_path = injection_name+ "_full.root"
    #print("full_root_path ",full_root_path)
    # Check if "_full.root" files exist
    if os.path.exists(full_root_path) or (os.path.exists(injection_full_list_path) and full_root_path in open(injection_full_list_path).read()):
        #print("The data files has been processed:")
        #print(full_root_path)
        if os.path.exists(injection_full_list_path):
            if full_root_path not in open(injection_full_list_path).read():
                with open(injection_full_list_path, "a") as file:
                    file.write(full_root_path + "\n")
        return True
    else:
        return False
    
                                
def create_and_fill_TH2F(injection_list_path, injection_name, file_path, frameID, frameID_start, frameID_stop):
    if not check_root_data_files_exist(injection_list_path, injection_name):
        if  check_bin_data_files_exist(file_path):
            # Call the function
            fre, time, amp = adjust_and_load_frames(file_path, frameID, frameID_start, frameID_stop)
            # Calculate the bins number
            nSpecBins, nTimeBins = len(fre), len(time)
            # Reshape the 1D array to a 2D array with dimensions (nSpecBins, nTimeBins)
            #amp = amp.reshape((nTimeBins, nSpecBins))
            
            # Calculate the boundaries of the bins
            freMin, freMax = min(fre), max(fre)
            timeMin, timeMax = min(time), max(time)
            
            # Assuming a uniform distribution, calculate the bin widths
            freBinWidth = (freMax - freMin) / (nSpecBins - 1)
            #timeBinWidth = (timeMax - timeMin) / (nTimeBins - 1)
            #timeBinWidth = (time[1] - time[0])
            timeBinWidth = 1./20000000*262144
            
            # Create the TH2F histogram
            last_part = injection_name.split('/')[-1]
            hist_name = "FFT_injection"
            hist_title = last_part +";Frequency (kHz) + 245 MHz;Time (s)"
            x1 = (freMin  -  freBinWidth / 2)/1000.
            x2 = (freMax  +  freBinWidth / 2)/1000.
            
            y1 = (frameID_start - frameID - 0.5)*timeBinWidth
            y2 = (frameID_stop  - frameID + 0.5)*timeBinWidth
            
            fHist0_ = TH2F(hist_name, hist_title, nSpecBins, x1, x2, nTimeBins, y1, y2)
            # Assuming 'amp' is a 2D array where amp[i][j] corresponds to the amplitude
            # for the ith frequency and jth time
            #with open("data.txt", "w") as file:  # 打开文件data.txt用于写入，如果文件不存在则创建
            for j in range(len(time)):               
                for i in range(len(fre)):
                    y = fHist0_.GetYaxis().GetBinCenter(j)
                    fHist0_.Fill(fre[i]/1000., y, amp[j][i])
            #print("chenrj injection_list_path = ",injection_list_path)        
            name = injection_name + "_full.root"
            f = TFile(name,"recreate")
            fHist0_.Write()
            f.Close()
            injection_full_list_path = injection_list_path.replace("_list.txt","_full_list.txt")
            # Append 'name' to the end of the injection_full_list_path file
            with open(injection_full_list_path, "a") as file:
                file.write(name + "\n")
        else:
            print("One or more required data files are missing.")
            
def process_injection(injection_list_path, injection_data, base_path, frameID_offset, frameID_range):
    """
    Process a single injection. This function is intended to be run in a separate process.
    
    Parameters:
    - injection_data: A tuple containing (injection_name, file_path, frameID).
    - base_path: The base path to replace the relative data path prefix.
    """
    injection_name, file_path, frameID = injection_data
    
    print(f"####### Processing {injection_name}, file_path = {file_path}, frameID = {frameID}")
    injection_name = injection_name.replace("../data", base_path)
    injection_name = injection_name.replace(".root", "")
    
    # Call the function to create and fill the histogram
    create_and_fill_TH2F(injection_list_path, injection_name, file_path, frameID, frameID - frameID_offset, frameID - frameID_offset+ frameID_range)

def convert_to_time_offset(date_start):
    # Convert the date_start string to a datetime object
    dt_start = datetime.strptime(date_start, "%Y-%m-%d %H:%M:%S")
    
    # Create a TDatime object using the datetime values
    tdatime_start = TDatime(dt_start.year, dt_start.month, dt_start.day, dt_start.hour, dt_start.minute, dt_start.second)
    
    # Get the time offset in seconds
    time_offset = tdatime_start.Convert() - 25*365*24*3600 -6*24*3600
    
    return time_offset

            
def combine_tdms(file_dir, file_start, file_stop, average_number, fre_min, fre_range,  date_start):
    # Initialize lists to accumulate the filtered data
    accumulated_fre = []
    accumulated_time = []
    accumulated_amp = []
    
    time_offset = convert_to_time_offset(date_start)
    
    # Process each file in the list
    i = 0
    fre_size = None
    fre      = None
    time     = None
    amp      = None
    mask = None
    fre_filtered = None
    for file_index in range(file_start, file_stop + 1):
        print("Processing file: ", file_index)
        file_name_prefix = f"{file_dir}/{file_index:07d}.iq.tdms"
        bin_time_path    = f"{file_name_prefix}.bin_time"
        bin_fre_path     = f"{file_name_prefix}.bin_fre"
        bin_amp_path     = f"{file_name_prefix}.bin_amp"
        try:
            if i ==0:
                fre  = np.fromfile(bin_fre_path, dtype=np.float64)
                time = np.fromfile(bin_time_path, dtype=np.float32)
                
            amp  = np.fromfile(bin_amp_path, dtype=np.float32)
            
            # Ensure amp is reshaped correctly
            if amp.size % fre.size != 0 or amp.size != fre.size * time.size:
                raise ValueError("The size of amp is not correct, indicating data loss. Filling amp with 1e-2.")
            
            num_time_points = amp.size // fre.size
            amp = amp.reshape(num_time_points, fre.size)
            
            #except FileNotFoundError:
        except (FileNotFoundError, ValueError) as e:
            print(f"File {file_name_prefix} not found. Setting amp to 0.")
            accumulated_time.append(time + time_offset + file_index * max(time))
            accumulated_amp.append(np.zeros((time.size, fre_size)) + 1e-2)
            continue

        # Filter the fre and amp arrays based on the given range
        if i ==0:
            mask = (fre >= fre_min) & (fre <= fre_min + fre_range)
            fre_filtered = fre[mask]
            
        amp_filtered = amp[:, mask]
            
        if fre_size is None:
            fre_size = len(fre_filtered)
        elif fre_size != len(fre_filtered):
            raise ValueError("Frequency size mismatch. All files should have the same frequency dimensions after filtering.")
        
        # Accumulate the filtered data
        if (i == 0):
            accumulated_fre.append(fre_filtered)
            
        accumulated_time.append(time + time_offset + file_index*max(time))
        accumulated_amp.append(amp_filtered)
        i = i + 1
        
    # Concatenate all accumulated data
    total_fre  = np.concatenate(accumulated_fre)
    total_time = np.concatenate(accumulated_time)
    total_amp  = np.concatenate(accumulated_amp, axis=0)
    
    # Combine rows along the y-axis (time axis)
    num_combined_time_points = total_time.size // average_number
    combined_amp = np.zeros((num_combined_time_points, total_fre.size))
    
    for i in range(num_combined_time_points):
        start_index = i * average_number
        end_index = start_index + average_number
        combined_amp[i, :] = np.mean(total_amp[start_index:end_index, :], axis=0)
        
    combined_time = np.mean(total_time.reshape(-1, average_number), axis=1)
    
    bin_npz_path  = f"{file_dir}/filtered_spectrum_{file_start}_{file_stop}.npz"
    print(f"saving data into {bin_npz_path} file.")
    np.savez_compressed(bin_npz_path, arr_0=total_fre, arr_1=combined_time, arr_2=combined_amp)
    
    # Create a 2D histogram
    h_spec_filtered = TH2F("h_spec_filtered", "Filtered Spectrum",
                           num_combined_time_points, min(combined_time), max(combined_time),
                           len(total_fre), min(total_fre), max(total_fre))
    
    # Fill the histogram with filtered data
    for i in range(len(total_fre)):
        for j in range(num_combined_time_points):
            h_spec_filtered.Fill(combined_time[j], total_fre[i], combined_amp[j, i])
                    
    # Set time display format for X-axis
    h_spec_filtered.GetXaxis().SetTimeDisplay(1)
    h_spec_filtered.GetXaxis().SetTimeFormat("%Y-%m-%d %H:%M:%S")
    h_spec_filtered.GetXaxis().SetLabelSize(0.05)
    h_spec_filtered.GetXaxis().SetTitleSize(0.05)
    h_spec_filtered.GetXaxis().SetTitle("Date")
    h_spec_filtered.GetXaxis().SetNdivisions(505)
    h_spec_filtered.GetYaxis().SetTitle("Frequency [Hz]")
    h_spec_filtered.GetYaxis().SetLabelSize(0.05)
    h_spec_filtered.GetYaxis().SetTitleSize(0.05)
                            
    # Save the histogram to a ROOT file
    output_root_file = f"{file_dir}/filtered_spectrum_{file_start}_{file_stop}.root"
    root_file = TFile(output_root_file, "RECREATE")
    h_spec_filtered.Write()
    root_file.Close()
    print(f"Histogram saved in {output_root_file}.")

    # Create a canvas and draw the histogram
    c1 = TCanvas("c1", "Filtered Spectrum", 1200, 600)
    c1.Divide(1, 2)
    # Adjust margins
    c1.cd(1).SetTopMargin(0.0)
    c1.cd(1).SetRightMargin(0.0)
    c1.cd(2).SetTopMargin(0.0)
    c1.cd(2).SetRightMargin(0.0)
    
    # Draw the 2D histogram
    c1.cd(1)
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gPad.SetLogz()
    h_spec_filtered.Draw("COLZ")
    
    # Draw the projection on the x-axis
    c1.cd(2)
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)
    gPad.SetLogy()
    h_proj_time = h_spec_filtered.ProjectionX()
    h_proj_time.SetTitle("Projection on Time Axis;Date;Amplitude Sum")
    h_proj_time.GetXaxis().SetTimeDisplay(1)
    h_proj_time.GetXaxis().SetTimeFormat("%Y-%m-%d %H:%M:%S")
    h_proj_time.GetXaxis().SetLabelSize(0.05)
    h_proj_time.GetXaxis().SetTitleSize(0.05)
    h_proj_time.GetXaxis().SetNdivisions(505)
    h_proj_time.GetYaxis().SetLabelSize(0.05)
    h_proj_time.GetYaxis().SetTitleSize(0.05)
    h_proj_time.GetYaxis().SetNdivisions(505)
    h_proj_time.Draw()
    # Show the canvas
    c1.Update()
    c1.Draw()
    gApplication.Run(True)  # Keep the canvas open for further editing
    
class CombineInjectionGUI(QWidget):
    def __init__(self):
        super().__init__()
        
        self.config_file = "parameters.toml"
        self.initUI()
        self.load_parameters()
        
    def initUI(self):
        self.setWindowTitle('Combine Injection GUI')
        self.setGeometry(100, 100, 600, 400)
        
        self.layout = QVBoxLayout()
    
        # Mode selection
        self.mode_label = QLabel("Mode:")
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["accumulate", "get_full_spec", "convert_npz", "combine_tdms"])
        self.mode_combo.currentIndexChanged.connect(self.update_ui_for_mode)
        self.layout.addWidget(self.mode_label)
        self.layout.addWidget(self.mode_combo)
        
        # Parameters layout
        self.param_layout = QFormLayout()
        self.layout.addLayout(self.param_layout)
        
        self.start_button = QPushButton("Start", self)
        self.start_button.clicked.connect(self.start_program)
        self.layout.addWidget(self.start_button)
        
        self.setLayout(self.layout)
        self.update_ui_for_mode()
        
    def update_ui_for_mode(self):
        print("Updating UI for mode:", self.mode_combo.currentText())
        # Clear previous parameters
        count = self.param_layout.count()
        print(f"Parameter layout count before clearing: {count}")
        while count > 0:
            child = self.param_layout.takeAt(0)
            if child is not None:
                widget = child.widget()
                if widget is not None:
                    print("Removing widget:", widget.objectName())
                    widget.setParent(None)
                    widget.deleteLater()
            count -= 1
                
        count_after = self.param_layout.count()
        print(f"Parameter layout count after clearing: {count_after}")
        
        mode = self.mode_combo.currentText()
        self.param_entries = []
        
        if mode == "accumulate":
            self.create_param_entry("path_to_files")
            self.create_param_entry("prefix")
            self.create_param_entry("extension")
            self.create_param_entry("start_time")
            self.create_param_entry("end_time")
            self.create_param_entry("dtime")
        elif mode == "get_full_spec":
            self.create_param_entry("injection_list_path")
            self.create_param_entry("base_path")
            self.create_param_entry("max_workers")
            self.create_param_entry("frameID_offset")
            self.create_param_entry("frameID_range")
        elif mode == "convert_npz":
            self.create_param_entry("bin_fre_path")
            self.create_param_entry("bin_time_path")
            self.create_param_entry("bin_amp_path")
        elif mode == "combine_tdms":
            self.create_param_entry("file_dir")
            self.create_param_entry("file_start")
            self.create_param_entry("file_stop")
            self.create_param_entry("average_number")
            self.create_param_entry("fre_min")
            self.create_param_entry("fre_range")
            self.create_param_entry("date_start")
        
    def create_param_entry(self, param_name):
        label = QLabel(f"{param_name}:")
        entry = QLineEdit(self)
        entry.setObjectName(param_name)
        self.param_entries.append((param_name, entry))
        self.param_layout.addRow(label, entry)
        print(f"Created entry for {param_name}")
        
    def load_parameters(self):
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as file:
                config = toml.load(file)
                mode = config.get("mode", "accumulate")
                self.mode_combo.setCurrentText(mode)
                params = config.get("parameters", {})
                for param_name, entry in self.param_entries:
                    entry.setText(params.get(param_name, ""))
                        
    def save_parameters(self):
        params = {param_name: entry.text() for param_name, entry in self.param_entries}
        config = {"mode": self.mode_combo.currentText(), "parameters": params}
        with open(self.config_file, 'w') as file:
            toml.dump(config, file)
                
    def start_program(self):
        self.save_parameters()
        params = {param_name: entry.text() for param_name, entry in self.param_entries}
        mode = self.mode_combo.currentText()
        
        if mode == "accumulate":
            path_to_files = params["path_to_files"]
            prefix = params["prefix"]
            extension = params["extension"]
            start_time = params["start_time"]
            end_time = params["end_time"]
            dtime = float(params["dtime"])
            filtered_files = filter_files_by_time(path_to_files, prefix, extension, start_time, end_time, "%Y.%m.%d.%H.%M.%S.%f")
            read_fft_injection_from_files(filtered_files, extension, dtime)
        elif mode == "get_full_spec":
            injection_list_path = params["injection_list_path"]
            base_path = params["base_path"]
            max_workers = int(params["max_workers"])
            frameID_offset = int(params["frameID_offset"])
            frameID_range = int(params["frameID_range"])
            injections = read_injection_list(injection_list_path)
            process_injections_concurrently(injection_list_path, injections, base_path, process_injection, max_workers, frameID_offset, frameID_range)
        elif mode == "convert_npz":
            bin_fre_path = params["bin_fre_path"]
            bin_time_path = params["bin_time_path"]
            bin_amp_path = params["bin_amp_path"]
            start_time = time.time()
            bin_npz_path = bin_fre_path.replace(".bin_fre", ".npz")
            fre = np.fromfile(bin_fre_path, dtype=np.float64)
            time_data = np.fromfile(bin_time_path, dtype=np.float32)
            amp = np.fromfile(bin_amp_path, dtype=np.float32)
            print("saving data into npz file.")
            np.savez_compressed(bin_npz_path, fre=fre, time=time_data, amp=amp)
            print(f"Data saved to {bin_npz_path}")
            end_time = time.time()
            used_time = end_time - start_time
            print(f"Time used: {used_time} seconds")
        elif mode == "combine_tdms":
            file_dir = params["file_dir"]
            file_start = int(params["file_start"])
            file_stop = int(params["file_stop"])
            average_number = int(params["average_number"])
            fre_min = float(params["fre_min"])
            fre_range = float(params["fre_range"])
            date_start = params["date_start"]
            combine_tdms(file_dir, file_start, file_stop, average_number, fre_min, fre_range, date_start)
            
def main():
    if len(sys.argv) > 1 and sys.argv[1].endswith(".toml"):
        config_file = sys.argv[1]
        app = QApplication(sys.argv)
        gui = CombineInjectionGUI()
        gui.config_file = config_file
        gui.load_parameters()
        
        #a1gui.start_program()
        gui.show()
        sys.exit(app.exec_())
    else:
        app = QApplication(sys.argv)
        gui = CombineInjectionGUI()
        gui.show()
        sys.exit(app.exec_())
