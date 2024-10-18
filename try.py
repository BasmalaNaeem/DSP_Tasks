import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
import numpy as np
import matplotlib.pyplot as plt
def ReadSignalFile(file_name):
    expected_indices = []
    expected_samples = []
    with open(file_name, 'r') as f:
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        while line:
            # Process line
            L = line.strip()
            if len(L.split(' ')) == 2:
                L = line.split(' ')
                V1 = int(L[0])
                V2 = float(L[1])
                expected_indices.append(V1)
                expected_samples.append(V2)
                line = f.readline()
            else:
                break
    return expected_indices, expected_samples



def plot_signal(indices, samples, title="Signal", color="b"):
    plt.figure()
    plt.plot(indices, samples, color=color)
    plt.title(title)
    plt.xlabel("Index")
    plt.ylabel("Sample Value")
    plt.grid(True)
    plt.show()


def save_signal_to_file(signal, file_path):
    with open(file_path, 'w') as f:
        f.write(f"{len(signal)}\n")
        
        # Write each index and sample
        for index, sample in signal:
            # Ensure the sample is written is int
            if isinstance(sample, float) and sample.is_integer():
                f.write(f"{index} {int(sample)}\n")
            else:
                f.write(f"{index} {sample}\n")



root = tk.Tk()
root.title("Signal Processing")


def AddSignalSamplesAreEqual(userFirstSignal,userSecondSignal,Your_indices,Your_samples):
    if(userFirstSignal=='Signal1.txt' and userSecondSignal=='Signal2.txt'):
        file_name="add.txt"  # write here the path of the add output file
    expected_indices,expected_samples=ReadSignalFile(file_name)          
    if (len(expected_samples)!=len(Your_samples)) and (len(expected_indices)!=len(Your_indices)):
        print("Addition Test case failed, your signal have different length from the expected one")
        return
    for i in range(len(Your_indices)):
        if(Your_indices[i]!=expected_indices[i]):
            print("Addition Test case failed, your signal have different indicies from the expected one") 
            return
    for i in range(len(expected_samples)):
        if abs(Your_samples[i] - expected_samples[i]) < 0.01:
            continue
        else:
            print("Addition Test case failed, your signal have different values from the expected one") 
            return
    print("Addition Test case passed successfully")

def load_and_plot_signal():
    file_name1 = "Signal1.txt"
    file_name2 = "Signal2.txt"
    
    indices1, samples1 = ReadSignalFile(file_name1)
    plot_signal(indices1, samples1, title="Signal 1")

    indices2, samples2 = ReadSignalFile(file_name2)
    plot_signal(indices2, samples2, title="Signal 2")

def add_signals():
    file_name1 = "Signal1.txt"
    file_name2 = "Signal2.txt"
    
    indices1, samples1 = ReadSignalFile(file_name1)
    indices2, samples2 = ReadSignalFile(file_name2)
    
    combined_indices = sorted(set(indices1) | set(indices2))
    combined_samples = []
    
    for index in combined_indices:
        value1 = samples1[indices1.index(index)] if index in indices1 else 0
        value2 = samples2[indices2.index(index)] if index in indices2 else 0
        combined_value = value1 + value2
        
        if combined_value.is_integer():
            combined_samples.append(int(combined_value))
        else:
            combined_samples.append(combined_value)

    plot_signal(combined_indices, combined_samples, title="Added Signal", color="g")
    
    AddSignalSamplesAreEqual(file_name1, file_name2, combined_indices, combined_samples)

    save_file_path = "added_signal.txt"
    save_signal_to_file(list(zip(combined_indices, combined_samples)), save_file_path)
    messagebox.showinfo("Success", "Added signal saved successfully!")


def SubSignalSamplesAreEqual(userFirstSignal,userSecondSignal,Your_indices,Your_samples):
    if(userFirstSignal=='Signal1.txt' and userSecondSignal=='Signal2.txt'):
        file_name="subtract.txt" # write here the path of the subtract output file
        
    expected_indices,expected_samples=ReadSignalFile(file_name)   
    
    if (len(expected_samples)!=len(Your_samples)) and (len(expected_indices)!=len(Your_indices)):
        print("Subtraction Test case failed, your signal have different length from the expected one")
        return
    for i in range(len(Your_indices)):
        if(Your_indices[i]!=expected_indices[i]):
            print("Subtraction Test case failed, your signal have different indicies from the expected one") 
            return
    for i in range(len(expected_samples)):
        if abs(Your_samples[i] - expected_samples[i]) < 0.01:
            continue
        else:
            print("Subtraction Test case failed, your signal have different values from the expected one") 
            return
    print("Subtraction Test case passed successfully")

def subtract_signals():
    file_name1 = "Signal1.txt"
    file_name2 = "Signal2.txt"
    
    indices1, samples1 = ReadSignalFile(file_name1)
    indices2, samples2 = ReadSignalFile(file_name2)
    
    combined_indices = sorted(set(indices1) | set(indices2))
    combined_samples = []
    
    for index in combined_indices:
        value1 = samples1[indices1.index(index)] if index in indices1 else 0
        value2 = samples2[indices2.index(index)] if index in indices2 else 0
        combined_value = value1 - value2
        
        if combined_value.is_integer():
            combined_samples.append(int(combined_value))
        else:
            combined_samples.append(combined_value)

    plot_signal(combined_indices, combined_samples, title="Subtracted Signal", color="r")
    SubSignalSamplesAreEqual(file_name1, file_name2,combined_indices, combined_samples) 

    save_file_path = "subtracted_signal.txt"  
    save_signal_to_file(list(zip(combined_indices, combined_samples)), save_file_path)
    messagebox.showinfo("Success", "Subtracted signal saved successfully!")

def ShiftSignalByConst(Shift_value,Your_indices,Your_samples):
    if(Shift_value==3):  #x(n+k)
        file_name="advance3.txt" # write here the path of delay3 output file
    elif(Shift_value==-3): #x(n-k)
        file_name="delay3.txt" # write here the path of advance3 output file
        
    expected_indices,expected_samples=ReadSignalFile(file_name)      
    if (len(expected_samples)!=len(Your_samples)) and (len(expected_indices)!=len(Your_indices)):
        print("Shift by "+str(Shift_value)+" Test case failed, your signal have different length from the expected one")
        return
    for i in range(len(Your_indices)):
        if(Your_indices[i]!=expected_indices[i]):
            print("Shift by "+str(Shift_value)+" Test case failed, your signal have different indicies from the expected one") 
            return
    for i in range(len(expected_samples)):
        if abs(Your_samples[i] - expected_samples[i]) < 0.01:
            continue
        else:
            print("Shift by "+str(Shift_value)+" Test case failed, your signal have different values from the expected one") 
            return
    print("Shift by "+str(Shift_value)+" Test case passed successfully")

def shift_signal():
    file_name = "Signal1.txt"
    
    original_indices, original_samples = ReadSignalFile(file_name)
                                                                    # -1   x(n+1) shfit by 1 advance - x axis
    shift_value = simpledialog.askinteger("Input", "Enter shift constant (positive for advancing, negative for delaying):")
    
    if shift_value is not None:
        shifted_indices = [index - shift_value for index in original_indices]
        
        shifted_samples = original_samples
        
        plot_signal(shifted_indices, shifted_samples, title=f"Signal Shifted by {shift_value}", color="b")
        
        ShiftSignalByConst(shift_value, shifted_indices, shifted_samples) 
        
        save_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if save_file_path:
            shifted_signal = list(zip(shifted_indices, shifted_samples))
            save_signal_to_file(shifted_signal, save_file_path)
            messagebox.showinfo("Success", "Shifted signal saved successfully!")

def Folding(Your_indices,Your_samples):
    file_name = "folding.txt"  # write here the path of the folding output file
    expected_indices,expected_samples=ReadSignalFile(file_name)      
    if (len(expected_samples)!=len(Your_samples)) and (len(expected_indices)!=len(Your_indices)):
        print("Folding Test case failed, your signal have different length from the expected one")
        return
    for i in range(len(Your_indices)):
        if(Your_indices[i]!=expected_indices[i]):
            print("Folding Test case failed, your signal have different indicies from the expected one") 
            return
    for i in range(len(expected_samples)):
        if abs(Your_samples[i] - expected_samples[i]) < 0.01:
            continue
        else:
            print("Folding Test case failed, your signal have different values from the expected one") 
            return
    print("Folding Test case passed successfully")

def fold_signal():
    file_name = "Signal1.txt"
    
    original_indices, original_samples = ReadSignalFile(file_name)
    
    folded_signal = [(index * -1, value) for index, value in zip(original_indices, original_samples)]
    
    folded_signal.sort(key=lambda x: x[0])
    
    folded_indices = [index for index, _ in folded_signal]
    folded_samples = [value for _, value in folded_signal]
    
    plot_signal(folded_indices, folded_samples, title="Folded/Reversed Signal x(-n)", color="r")
    
    Folding(folded_indices, folded_samples)
    
    save_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
    if save_file_path:
        save_signal_to_file(folded_signal, save_file_path)
        messagebox.showinfo("Success", "Folded/Reversed signal saved successfully!")


def MultiplySignalByConst(User_Const,Your_indices,Your_samples):
    if(User_Const==5):
        file_name="mul5.txt"  # write here the path of the mul5 output file
        
    expected_indices,expected_samples=ReadSignalFile(file_name)      
    if (len(expected_samples)!=len(Your_samples)) and (len(expected_indices)!=len(Your_indices)):
        print("Multiply by "+str(User_Const)+ " Test case failed, your signal have different length from the expected one")
        return
    for i in range(len(Your_indices)):
        if(Your_indices[i]!=expected_indices[i]):
            print("Multiply by "+str(User_Const)+" Test case failed, your signal have different indicies from the expected one") 
            return
    for i in range(len(expected_samples)):
        if abs(Your_samples[i] - expected_samples[i]) < 0.01:
            continue
        else:
            print("Multiply by "+str(User_Const)+" Test case failed, your signal have different values from the expected one") 
            return
    print("Multiply by "+str(User_Const)+" Test case passed successfully")

def multiply_signal_with_constant():
    file_name = "Signal1.txt"
    
    original_indices, original_samples = ReadSignalFile(file_name)
    
    constant = simpledialog.askfloat("Input", "Enter multiplication constant:", minvalue=-10.0, maxvalue=10.0)
    
    if constant is not None:
        multiplied_signal = []
        
        for index, value in zip(original_indices, original_samples):
            new_value = value * constant 
            multiplied_signal.append((index, new_value)) 
            
        multiplied_indices = [sample[0] for sample in multiplied_signal]
        multiplied_samples = [sample[1] for sample in multiplied_signal]
        
        plot_signal(multiplied_indices, multiplied_samples, title=f"Signal Values Multiplied by {constant}", color="g")
        
        MultiplySignalByConst(constant, multiplied_indices, multiplied_samples)
        
        save_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        if save_file_path:
            save_signal_to_file(multiplied_signal, save_file_path)
            messagebox.showinfo("Success", "Multiplied signal saved successfully!")




button_style = {
    'bg': '#4DD0E1',        
    'fg': 'white',          
    'font': ('Helvetica', 12, 'bold'),  
    'width': 25,            
    'pady': 10,           
    'borderwidth': 3 
}

# load_button = tk.Button(root, text="Load Signal", command=loadandplotone, **button_style)
# load_button.pack(pady=5)
# advance_button = tk.Button(root, text="Advance Signal", command=advance_signal, **button_style)
# advance_button.pack(pady=5)

# delay_button = tk.Button(root, text="Delay Signal", command=delay_signal, **button_style)
# delay_button.pack(pady=5)

load_button_2 = tk.Button(root, text="Load Signals", command=load_and_plot_signal, **button_style)
load_button_2.pack(pady=10)

add_button = tk.Button(root, text="Add Signals", command=add_signals, **button_style)
add_button.pack(pady=5)

subtract_button = tk.Button(root, text="Subtract Signals", command=subtract_signals, **button_style)
subtract_button.pack(pady=5)

advance_button = tk.Button(root, text="Shift Signal", command=shift_signal, **button_style)
advance_button.pack(pady=5)

fold_button = tk.Button(root, text="Fold Signal", command=fold_signal, **button_style)
fold_button.pack(pady=5)

multiply_button = tk.Button(root, text="Multiply Signal by Constant", command=multiply_signal_with_constant, **button_style)
multiply_button.pack(pady=5)

root.mainloop()
