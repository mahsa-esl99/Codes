# DSP Lab Final Project
# N18123130
# N14334042

from math import pi, sin
import pyaudio
import struct
import tkinter as Tk
from PIL import ImageTk, Image
import scipy
from tkinter import messagebox
import wave
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import  butter, lfilter
from skimage.restoration import denoise_wavelet

global CONTINUE
CONTINUE = True
global Denoising
Denoising = False


def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

# Function for Quit button
def fun_quit():
    global CONTINUE
    print('Good bye')
    CONTINUE = False
    root.destroy()


# Function for de-noising the input signal
def denoise_func():
    global Denoising
    print('Noise Cancellation is in Progress ...')
    Denoising = True


# Define Tkinter root
root = Tk.Tk()
root.title('ANC GUI')
root.geometry("400x350")
frame = Tk.Frame(root, width=200, height=100)
frame.pack()
frame.place(anchor='center', relx=0.5, rely=0.4)

img = ImageTk.PhotoImage(Image.open("logo.jpg"))
label = Tk.Label(frame, image=img)
label.pack()

# Define widgets
#messagebox.showinfo('Info', 'Welcome to the demo of our project!')
B_noise = Tk.Button(root, text='Start Noise Cancellation', command=denoise_func, height=3, width=30)
B_quit = Tk.Button(root, text='Quit', command=fun_quit, height=3, width=5)

# Place widgets
B_quit.place(x=25, y=50)
B_noise.place(x=90, y=270)

WIDTH = 2           # bytes per sample
CHANNELS = 1        # mono
RATE = 8000     # frames per second
RECORD_SECONDS = 0.1     # Duration in seconds
BLOCKLEN = 128
MAXVALUE = 2**15-1  # Maximum allowed output signal value (because WIDTH = 2)
num_blocks = int(RATE / BLOCKLEN * RECORD_SECONDS)
order = 7
cutoff = 3500  # desired cutoff frequency of the filter, Hz

# Get the filter coefficients so we can check its frequency response.
b, a = butter_lowpass(cutoff, RATE, order)

p = pyaudio.PyAudio()

output_wavefile = 'Output_file.wav'
wf = wave.open(output_wavefile, 'w')      # Wave file for saving the siren sound
wf.setnchannels(1)                        # One channel (mono)
wf.setsampwidth(2)                        # Two bytes per sample
wf.setframerate(RATE)                       # Samples per second

# Open audio stream
stream = p.open(
    format=p.get_format_from_width(WIDTH),
    channels=CHANNELS,
    rate=RATE,
    input=True,
    output=True)


#plt.ion()           # Turn on interactive mode so plot gets updated
DBscale = False
#DBscale = True

# Initialize plot window:
#plt.figure(1)
# if DBscale:
#     plt.ylim(0, 150)
# else:
#     plt.ylim(0, 5*RATE)

# Frequency axis (Hz)
# plt.xlim(0, 0.5*RATE)         # set x-axis limits
# plt.xlim(0, 2000)         # set x-axis limits
# plt.xlabel('Frequency (Hz)')
f = RATE/BLOCKLEN * np.arange(0, BLOCKLEN)

# line1, = plt.plot([], [], color='blue')  # Create empty line
# line1.set_xdata(f)                         # x-data of plot (frequency)
# line2, = plt.plot([], [], color='red')  # Create empty line
# line2.set_xdata(f)                         # x-data of plot (frequency)
# line1.set_label('Input')  # Set parameters of line1
# line2.set_label('Output')  # Set parameters of line2
# plt.legend()      # Create legend

print('* Start')
while CONTINUE:
    root.update()
    for i in range(0, num_blocks):
        
        input_bytes = stream.read(BLOCKLEN, exception_on_overflow=False)  # BLOCKLEN = number of frames read
        
        # Convert binary data to tuple of numbers
        input_tuple = struct.unpack('h' * BLOCKLEN, input_bytes)
        X = np.fft.fft(input_tuple)
        # f, t, X = scipy.signal.stft(input_tuple, RATE)
        # X = scipy.fft.fftshift(input_tuple)
        S = np.ndarray(shape=(len(X)), dtype=complex)
        Y = np.ndarray(shape=(len(X)), dtype=complex)
        # if DBscale:
        #     line1.set_ydata(20 * np.log10(np.abs(X)))
        # else:
        #     line1.set_ydata(np.abs(X))
        if Denoising:
            var_noise = 0.15*10**8
            var_temp = np.mean(np.abs(X**2)) - var_noise
            
            if var_temp >= 0:
                var_est = np.sqrt(var_temp)
            else:
                var_est = 10
                
            T = np.sqrt(2)*var_noise/var_est
            for k in range(len(X)):
                if (X[k]) < -T:
                    S[k] = X[k] + T
                elif (X[k]) >= -T and (X[k]) <= T:
                    S[k] = 0
                elif (X[k]) > T:
                    S[k] = X[k] - T
            s = np.fft.ifft(S)       
            s_wavelet = denoise_wavelet(np.real(s), method='BayesShrink', mode='soft', wavelet_levels=1, wavelet='db3', rescale_sigma='True')
            # threshold_value = 0.06*np.max(abs(X.all()))
            # for k in range(len(X)):
            #     print(X[k])
            #     if np.abs(X[k]) <= threshold_value:
            #         S[k] = 0 
            #     elif np.abs(X[k]) > threshold_value:
            #         S[k] = X[k]
        else:
            S = X
            # s_wavelet = np.fft.ifft(S)
            s = np.fft.ifft(S)
            

        
        # _, s = scipy.signal.istft(S, RATE)
        # s = scipy.fft.ifft(S)
        
       # convert to integer
        # output_block = s_wavelet.astype(int)
        output_block = s.astype(int)
        #line2.set_ydata(np.abs(S))
        #plt.pause(0.0001)
        # Convert output value to binary data
        output_bytes = struct.pack('h' * BLOCKLEN, *output_block)

        # Write binary data to audio output stream
        stream.write(output_bytes)
        wf.writeframes(output_bytes)


wf.close()
stream.stop_stream()
stream.close()
p.terminate()

plt.ioff()           # Turn off interactive mode
plt.show()           # Keep plot showing at end of program
plt.close()
print('* Finished')