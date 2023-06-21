#! /usr/bin/env python


### Welcome to BracketMaker, a program for generating chemical protein synthesis (CPS) bracket figures and optimizing CPS strategies.
### Created by Judah Evangelista, Kay Lab, University of Utah
### To run BracketMaker GUI, run this script, which contains all elements of the GUI and GUI-specific functions.
### The primary functions of BracketMaker (including figure-making, auto-bracket, and bracket scoring) are defined in a separate file "scripts/bracketmaker.py" and are imported upon launch.
### "scripts/bracketmaker.py" must be present or the GUI will not work from this script! (Ignore for the executable)

#The following changes the working directory to the folder in which the BracketMaker executable is stored in.
#If using the Python.exe executable (i.e., running from a script instead of using BracketMaker executable), this will be ignored.
if not 'python' in os.path.basename(sys.executable):
    os.chdir(os.path.dirname(sys.executable))
# Report the current working directory...this says where /bm_cache folder may be found and verifies the above line worked correctly.
print(f"Running in folder: {os.getcwd()}")

import tkinter as tk
from tkinter import scrolledtext,filedialog,simpledialog,colorchooser,messagebox
import re
import os
import sys
import shutil
import platform # for Windows/Mac differences
import importlib
import drawSvg as draw # Main module used by BracketMaker to generate SVG images
from PIL import Image,ImageTk # Fxns for displaying images in tkinter

# Load primary BracketMaker functions
from scripts import bracketmaker


### This script is organized as follows (use Find to quickly navigate)
### I...GLOBAL VARIABLES ON STARTUP
### II...GUI FUNCTION DEFS
### III...GUI WINDOW ASSEMBLY


### ----------------------------
### I...GLOBAL VARIABLES ON STARTUP
### ----------------------------

# These settings influence behavior of GUI
reload_bracketmaker_on_run = True # Re-import bracketmaker script each time the function is called in GUI; this allows changes to bracketmaker.py to be visualized in real time
open_gui = True # Open GUI when this script is run
auto_load_last = False # Automatically load the last bracket stored in the bm_cache folder

# DEFINE GUI WINDOW
window = tk.Tk()
window.title("BracketMaker")

# CREATE BRACKETMAKER WORKING FOLDER
os.makedirs('bm_cache', exist_ok=True)

# Load icon image; otherwise default tkinter icon will be used
try:
    makericon=tk.PhotoImage(file='bm_cache/makericon.png')
    window.iconphoto(False,makericon)
except:
    print(f'Missing logo file bm_cache/makericon.png')

# GUI VARIABLES
global_lastsvg = ''

img_preview_max_width = 800 #px
img_preview_max_height = 450 #px

# Settings for look/feel of GUI
normal_text_color = 'black'
inside_bracket_color = 'orange' # Inside parentheses color is linked to usr_highlight_color (same as AA highlight in Sequence portion of output SVG)
after_hyphen_color = 'grey'

# Toggles to show/hide parts of the GUI
show_preview_image=tk.IntVar(window) # This controls whether preview image should currently be visible. IntVar can be linked to a user setting 'checkbox'; this has been removed from GUI, and is default True. If this turns from True>False, the image frame will be cleared.
show_preview_image.set(True)
show_preview_image.trace('w',lambda *args: toggleShowPreview())

auto_preview_image=tk.IntVar(window) # This controls whether preview image automatically updates with any change in sequence or variables. IntVar can be linked to a user setting 'checkbox'; this has been removed from GUI, and is default True.
auto_preview_image.set(True)

show_image_settings = False # This controls whether the right-hand Figure Settings menu is expanded or collapsed upon initial launch.

show_autobracket=False # This controls whether the right-hand Auto-Bracket menu is expanded or collapsed upon initial launch.

# DEFINE VARIABLES FOR USER SETTINGS
# Initial segment colors
default_colors=(('#DF9282', '#926358'),('#EADA6C', '#A1964C'),('#AAE1DF', '#789D9C'),('#ADA9E4', '#77749C'),('#E7A4E0', '#9A748A'),('#E7A963', '#9F7646'),('#94D374', '#6C9A56'),('#AEC0E6', '#727E98'),('#C6A6DA', '#715D7F')) # Tuple containing fill,stroke; this will be parsed into separate lists of fill/stroke colors when swatches are loaded.

# Use IntVar() for checkboxes and numerical things, StringVar() for text
# For a tk variable to be responsive in real-time, "trace" it to a loadPreviewImage function; see the block at the beginning of GUI assembly
usr_align_segments_at_top = tk.IntVar()
usr_show_thioester_placeholder = tk.IntVar()
usr_thioester = tk.StringVar(window)
usr_show_segment_labels = tk.IntVar(window)
usr_show_special_aa_annotations = tk.IntVar(window)
usr_show_rxn_label_placeholder = tk.IntVar(window)
usr_rxn_label_placeholder_text = tk.StringVar(window)
usr_aa_start_number = tk.StringVar(window)
usr_fill_color_list=[]
usr_stroke_color_list=[]
usr_highlight_color=""
usr_poor_thioesters = tk.StringVar(window) # BracketMaker and AutoBracket use 'good_thioesters', which is the inverse of this
usr_highlight_poor_thioesters = tk.IntVar()
usr_highlight_non_cys_thiols = tk.IntVar()
usr_px_per_aa = tk.StringVar(window)
usr_segment_spacing = tk.StringVar(window)

# DEFAULT SETTINGS - Set initial values of all variables
usr_align_segments_at_top.set(True)
usr_show_thioester_placeholder.set(True)
usr_thioester.set("NHNH2")
usr_show_segment_labels.set(True)
usr_show_special_aa_annotations.set(True)
usr_colors = list(default_colors)*10
usr_highlight_color = 'firebrick'
usr_show_rxn_label_placeholder.set(True)
usr_rxn_label_placeholder_text.set('NCL')
usr_aa_start_number.set(str(1))
usr_poor_thioesters.set('I,K,L,T,V')
usr_highlight_poor_thioesters.set(True)
usr_highlight_non_cys_thiols.set(True)
usr_px_per_aa.set('3')
usr_segment_spacing.set('20')

# Default Save filetypes are .txt and .svg
usr_save_txt = tk.IntVar(window)
usr_save_txt.set(1)
usr_save_svg = tk.IntVar(window)
usr_save_svg.set(1)
usr_save_png = tk.IntVar(window)
usr_save_png.set(0)

# AutoBracket initial settings
usr_one_pot_desulfurization = tk.IntVar(window)
usr_one_pot_desulfurization.set(0)
usr_one_pot_internal_cys = tk.IntVar(window)
usr_one_pot_internal_cys.set(0)
usr_one_potnterm_cys = tk.IntVar(window)
usr_one_potnterm_cys.set(1)
usr_save_autobracket_file = tk.IntVar(window)
usr_save_autobracket_file.set(1)
usr_penalize_poor_thioesters = tk.IntVar(window)
usr_penalize_poor_thioesters.set(1)
usr_penalize_non_cys_thiols = tk.IntVar(window)
usr_penalize_non_cys_thiols.set(1)

# Initial bracket scores - start blank
var_maxpath = tk.StringVar(window)
var_maxpath.set('-')
var_avgpath = tk.StringVar(window)
var_avgpath.set('-')
var_steps = tk.StringVar(window)
var_steps.set('-')

# Initialize file navigation
current_file_index = 0
total_file_number = 1
var_currentfile = tk.StringVar(window)
var_currentfile.set(f"{current_file_index+1}")
var_totalfiles = tk.StringVar(window)
var_totalfiles.set(f'/{total_file_number}')
sequence_cache_list = []

# Default initial values for protein name and sequence fields
ent_protname_initial = txt_sequence_initial = ''
# If toggled, instead load these from cached file
if auto_load_last == True:
    tempfilepath = 'bm_cache/currentsequence.txt'
    if os.path.exists(tempfilepath) == True:
        with open(tempfilepath, 'r') as f:
            filecontents = f.readlines() # Gives list of lines
        # If .txt file contains header line beginning with >, add this to the Protein Name field and add the rest to text editor
        if filecontents[0][0] == '>':
            ent_protname_initial = filecontents[0][1:].strip()
            txt_sequence_initial = ''.join(filecontents[1:]).strip()
        # If starting > is not encountered, fields will be left blank


## -----------------------------
### II...GUI FUNCTION DEFS; listed alphabetically
## -----------------------------

def btnAligatorImport(): # Open an Aligator output .txt file, run Auto-Bracket on each line, and (optionally) save results to new file
    # If I make changes here, also change btnAutoBracket if relevant
    if reload_bracketmaker_on_run == True:
        importlib.reload(bracketmaker)
    # Check syntax of settings fields
    if checkAutoBracketSettings()==False:
        return
    numbertoshow = 1
    aligatorlinestoread=1000 # Note - add as user variable?
    cyspg1 = ent_internal_cys_pg.get().strip()
    cyspg2 = ent_nterm_cys_pg.get().strip()
    if cyspg2 == '':
        cyspg2 = cyspg1
    onepotdesulf=usr_one_pot_desulfurization.get()
    onepotinternalcys=usr_one_pot_internal_cys.get()
    onepotntermcys=usr_one_potnterm_cys.get()
    if usr_penalize_poor_thioesters.get()==True:
        goodthioesterslist=getGoodThioesters(usr_poor_thioesters.get())
    else:
        goodthioesterslist=getGoodThioesters('')
    if usr_penalize_non_cys_thiols.get()==True:
        noncysthiolpenalty=2
    else:
        noncysthiolpenalty=0
    saveautobracketfile=usr_save_autobracket_file.get()
    # PROMPT TO OPEN ALIGATOR FILE - File contents (readlines() result) are the input for the importAligator function
    filepath = filedialog.askopenfilename(
        initialdir="",
        title="Open Aligator file (<proteinname> All Strategies.txt)",
        filetypes=(("Text Files (.txt)", ".txt"),)
        )
    if not filepath:
        return
    # Prompt for Protein Name (numbers will be added)
    ProteinName = simpledialog.askstring(title='Protein Name',prompt='Protein Name:')
    if ProteinName == None: # i.e. if Cancel is pressed
        return
    # Run aligator brackets function, which returns a list of brackets and their scores
    report(f'Now opening Aligator file {filepath} (this may take a few seconds to minutes)')
    with open(filepath, 'r') as f:
        filecontents=f.readlines()
        BracketsAndScoresList = bracketmaker.getAligatorBrackets(filecontents, brackets_per_set=numbertoshow, lines_to_read=aligatorlinestoread, internal_cys_pg=cyspg1,n_term_cys_pg=cyspg2,good_thioesters=goodthioesterslist,non_cys_thiol_penalty=noncysthiolpenalty,one_pot_cys_deprotection=onepotinternalcys,one_pot_desulfurization=onepotdesulf,one_pot_nterm_cys_deprotection=onepotntermcys,sort_by_bracket_score=False,verbose=True) # This returns tuple pairs, ("Bracket",{ScoreDict})
        BracketList = [x[0] for x in BracketsAndScoresList]
    # Clear all currently loaded sequences
    clearSeqList()
    # Get current cache and add all strategies to it, then set current position to the first on the list
    global sequence_cache_list
    global current_file_index
    global total_file_number
    SeqCount = 0
    for Seq,Scores in BracketsAndScoresList:
        # Increment number for protein naming scheme
        SeqCount+=1
        numberextension = str(SeqCount).zfill(len(str(len(BracketsAndScoresList)))) # Funky, but this auto detects how many zeroes to add
        CurrentProtName = ProteinName+f' - Aligator Strategy {Scores["a"]} - Score {Scores["aligator"]}'
        # Create tuple with protein name and sequence and append to open sequences
        sequence_cache_list.append((CurrentProtName,Seq))
    # Set current position to 1 and set total; update all string variables and generate preview image
    current_file_index = 0
    total_file_number = len(sequence_cache_list)
    loadSequenceStrings(None,cache=False)
    # Report success
    report(f'Done - got brackets for {total_file_number} Aligator strategies')
    # PROMPT FOR FILE SAVE
    if saveautobracketfile==True:
        savefile = filedialog.asksaveasfile(
            mode='w',
            title = "Save bracket scores file (.tsv)",
            defaultextension='.tsv',
            filetypes=[('Tab-Separated Values (.tsv)','.tsv')]
            )
        if savefile=="": # If user hits Cancel, cancel this part only.
            return
        # Write scores to file
        with open(savefile.name,'w') as f:
            f.write('Full Bracket\tAligator Rank\tAligator Score\t# Segments\tMax Path\tAvg Path\tRxn Steps\n')
            for Bracket,score in BracketsAndScoresList:
                f.write(f'{Bracket}\t{score["a"]}\t{score["aligator"]}\t{score["segments"]}\t{score["maxpath"]}\t{score["avgpath"]}\t{score["steps"]}\n')
        # Report success
        report(f'Done - saved brackets for {total_file_number} Aligator strategies to {savefile.name}')
    # Toggle highlight of poor thioesters
    usr_highlight_poor_thioesters.set(usr_penalize_poor_thioesters.get())
    usr_highlight_non_cys_thiols.set(usr_penalize_non_cys_thiols.get())
    return

def btnAutoBracket(): # Run Auto-Bracket and (optionally) save results to file
    # If changes are made to this fxn, also change btnAligatorImport if relevant
    if reload_bracketmaker_on_run == True:
        importlib.reload(bracketmaker)
    global usr_one_pot_desulfurization
    global usr_one_pot_internal_cys
    global usr_one_potnterm_cys
    global usr_save_autobracket_file
    global usr_penalize_poor_thioesters
    global usr_poor_thioesters
    # Check text box format, and get input sequence to pass to fxn
    InputSequence = txt_sequence.get(1.0, tk.END).strip()
    if not "[]" in InputSequence:
        report('Format error! To run Auto-Bracket, insert empty [] tags between segments. White space is ignored.')
        return
    # Get protein name...but don't worry about it if there isn't one
    ProteinName = ent_protname.get()
    # Check syntax of settings boxes using separate fxn
    if checkAutoBracketSettings() == False:
        return
    # RETRIEVE AUTOBRACKET SETTINGS FROM GUI ELEMENTS - same as btnAligatorImport; if I make changes there also change them here
    numbertoshow = int(ent_number_of_brackets.get().strip())
    if numbertoshow==0:
        numbertoshow=1000000000000
    cyspg1 = ent_internal_cys_pg.get().strip()
    cyspg2 = ent_nterm_cys_pg.get().strip()
    if cyspg2 == '':
        cyspg2 = cyspg1
    onepotdesulf=usr_one_pot_desulfurization.get()
    onepotinternalcys=usr_one_pot_internal_cys.get()
    onepotntermcys=usr_one_potnterm_cys.get()
    if usr_penalize_poor_thioesters.get()==True:
        goodthioesterslist=getGoodThioesters(usr_poor_thioesters.get())
    else:
        goodthioesterslist=getGoodThioesters('')
    if usr_penalize_non_cys_thiols.get()==True:
        noncysthiolpenalty=2
    else:
        noncysthiolpenalty=0
    saveautobracketfile=usr_save_autobracket_file.get()
    # PROMPT FOR FILE SAVE
    if saveautobracketfile==True:
        savefile = filedialog.asksaveasfile(
            mode='w',
            title = "Save bracket scores file (.tsv)",
            defaultextension='.tsv',
            filetypes=[('Tab-Separated Values (.tsv)','.tsv')]
            )
        if savefile=="": # If user hits Cancel, cancel fxn
            return
    # Run auto-bracket function, which returns list of formatted brackets along with a dictionary of scores we can write to a file
    NumberOfSegments=len(re.split(r"\[\]",InputSequence))
    report(f'Now building brackets for {ProteinName} with {NumberOfSegments} segments (for proteins >12 segments this can take a long time)...')
    SortedList = bracketmaker.autoBracket(InputSequence, returnnumber=numbertoshow, internal_cys_pg=cyspg1, n_term_cys_pg=cyspg2, one_pot_desulfurization=onepotdesulf, one_pot_cys_deprotection=onepotinternalcys, one_pot_nterm_cys_deprotection=onepotntermcys, return_scores=False, good_thioesters=goodthioesterslist, non_cys_thiol_penalty=noncysthiolpenalty, verbose=True)
    # WRITE TO FILE IF TOGGLED
    if saveautobracketfile==True:
        with open(savefile.name,'w') as f:
            f.write('Full Bracket\tMax Path\tAvg Path\tRxn Steps\tSegments\n')
            for Bracket in SortedList:
                score = bracketmaker.scoreBracket(Bracket,good_thioesters=goodthioesterslist,non_cys_thiol_penalty=noncysthiolpenalty,typical_rxn_yield=0.5) # Dictionary containing all scores
                print("score dict below...")
                print(score)
                f.write(f'{Bracket}\t{score["maxpath"]}\t{score["avgpath"]}\t{score["steps"]}\t{score["segments"]}\n')
    # Clear all currently loaded sequences
    clearSeqList()
    # Get current cache and add all strategies to it, then set current position to the first on the list
    global sequence_cache_list
    global current_file_index
    global total_file_number
    SeqCount = 0
    if ProteinName == '':
        ProteinName = r'[Protein Name]'
    for Seq in SortedList:
        # Increment number for protein naming scheme
        SeqCount+=1
        numberextension = str(SeqCount).zfill(len(str(len(SortedList)))) # Funky, but this auto detects how many zeroes to add
        CurrentProtName = ProteinName+f' Strategy {numberextension}'
        # Create tuple with protein name and sequence and append to open sequences
        sequence_cache_list.append((CurrentProtName,Seq))
    # Set current position to 1 and set total; update all string variables and generate preview image
    current_file_index = 0
    total_file_number = len(sequence_cache_list)
    loadSequenceStrings(None,cache=False)
    # Report success
    report(f'Auto-Bracket complete! Showing {total_file_number} possible ligation orders.')
    if saveautobracketfile==True:
        report(f'Auto-Bracket complete! Saved {total_file_number} ligation orders to {savefile.name}')
    # Toggle highlight of poor thioesters
    usr_highlight_poor_thioesters.set(usr_penalize_poor_thioesters.get())
    usr_highlight_non_cys_thiols.set(usr_penalize_non_cys_thiols.get())
    return

def btnChangeSegmentColor(index): # Change color of the relevant segment using tkinter built-in color picker
    global usr_colors
    # Open color picker
    color_code = colorchooser.askcolor(title ="Highlight Color") # tuple ((r,g,b),'hex')
    # Set fill color and auto-generate a darker stroke color
    fillhex=color_code[1]
    strokehex=getDarkerColor(color_code[0])
    # Update the list of user colors at the current position
    usr_colors[index] = (fillhex,strokehex)
    # Reload preview image
    loadPreviewImage()

def btnLoadSettings():
    # Prompt to load a .bmsettings file
    filename = filedialog.askopenfilename(
        initialdir="C:/Users/MainFrame/Desktop/",
        title="Open BracketMaker Settings File",
        filetypes=(("Settings (.bmsettings)", ".bmsettings"),)
        )
    if not filename:
        return
    # bmsettings file is python code, simply setting each tkinter variable to its last recorded value.
    try:
        with open(filename,'r') as f:
            for line in f.readlines():
                # For some reason, colors won't load properly unless I do it like this. The other variables work fine. (Probably because colors list is not a tkinter dynamic variable like the others)
                if 'usr_colors' in line:
                    colorslist=eval(re.search(r'usr_colors\s*=\s*(.*)',line)[1])
                    global usr_colors
                    usr_colors=list(colorslist)
                else:
                    exec(line)
        report(f'Loaded settings file {filename}')
        # Re-load preview image, which should get correct color swatches as well
        loadPreviewImage()
        return
    except:
        report('Error loading settings file!')
        return

def btnNavClone(): # Navigation menu; copy current sequence and insert at the next position in the queue
    global current_file_index
    global total_file_number
    global sequence_cache_list
    # Cache the current sequence in place if it is not blank
    OldTxtSequence = txt_sequence.get(1.0, tk.END).strip()
    OldProtName = ent_protname.get().strip()
    if len(sequence_cache_list)==0: # Upon loading program cache will be blank
        sequence_cache_list.append('')
    sequence_cache_list[current_file_index] = (OldProtName,OldTxtSequence)
    # If this is not the last file, remember the rest of the queue
    rightsequences=[]
    if len(sequence_cache_list)>(current_file_index+1):
        rightsequences = sequence_cache_list[current_file_index+1:]
    # Record cache up to this point
    newcachelist = sequence_cache_list[0:current_file_index+1]
    # Clone sequence to new location in cache
    newcachelist.append((OldProtName,OldTxtSequence))
    # Add right sequences
    for x in rightsequences:
        newcachelist.append(x)
    # Set new cache and global variables
    sequence_cache_list = newcachelist[:]
    current_file_index = current_file_index+1
    total_file_number = len(sequence_cache_list)
    var_currentfile.set(f"{current_file_index+1}")
    var_totalfiles.set(f'/{total_file_number}')
    # Highlight and generate new figure
    syntaxHighlight()
    loadPreviewImage()
    report(f'Cloned strategy {OldProtName} to sequence {current_file_index+1}')

def btnNavClose(): # Navigation menu; remove the current sequence from the queue and load the next
    global current_file_index
    global total_file_number
    global sequence_cache_list
    OldProtName = ent_protname.get()
    # If this is the only sequence open, simply make it blank
    if total_file_number == 1:
        ent_protname.delete(0, tk.END)
        txt_sequence.delete(1.0, tk.END)
        return
    # If this is anything but the last sequence, delete and load the next one in the queue
    if current_file_index < (len(sequence_cache_list)-1):
        leftsequences=sequence_cache_list[0:current_file_index]
        rightsequences=sequence_cache_list[current_file_index+1:]
        sequence_cache_list = leftsequences+rightsequences
        #var_currentfile.set(f'{current_file_index+1}') # no need to change current file index
        total_file_number=len(sequence_cache_list)
        var_totalfiles.set(f'/{total_file_number}')
        # Run function to load the current sequence as if Next button was pressed; include syntax and preview image
        loadSequenceStrings(None,cache=False)
        report(f'Deleted sequence {OldProtName}')
        return
    # If this is the last sequence, delete and load the one behind it
    elif current_file_index == (len(sequence_cache_list)-1):
        sequence_cache_list = sequence_cache_list[0:(len(sequence_cache_list)-1)]
        current_file_index -= 1
        total_file_number = len(sequence_cache_list)
        var_currentfile.set(f'{current_file_index+1}')
        var_totalfiles.set(f'/{total_file_number}')
        loadSequenceStrings(None,cache=False)
        report(f'Deleted sequence {OldProtName}')
        return

def btnNavNextFile(): # Navigation menu; load the next sequence in queue
    global current_file_index
    global total_file_number
    old_file_index = current_file_index
    # If less than total file count, update number, update filenumber entry field and run its function
    if current_file_index+1 < total_file_number:
        current_file_index+=1
        var_currentfile.set(f'{current_file_index+1}')
        loadSequenceStrings(None,cacheindex=old_file_index)

def btnNavPrevFile(): # Navigation menu; load the previous sequence in queue
    global current_file_index
    # If above 1, update number, update filenumber entry field and run its function
    old_file_index = current_file_index
    if current_file_index > 0:
        current_file_index-=1
        var_currentfile.set(f'{current_file_index+1}')
        loadSequenceStrings(None,cacheindex=old_file_index)

def btnOpenList(): # Open list of brackets from a .tsv file (e.g., output of Auto-Bracket or Aligator Import)
    filepath = filedialog.askopenfilename(
        title="Open Bracket List (Tab-Separated Values, .tsv)",
        filetypes=(("Tab-Separated Values (.tsv)", ".tsv"),)
        )
    if not filepath:
        return
    global sequence_cache_list
    global current_file_index
    global total_file_number
    # Dialog popup; returns True to keep files, False to get rid of files
    # Ignore if only sequence is empty
    if len(sequence_cache_list)>0:
        keepfiles = messagebox.askyesnocancel("Question", "Do you want to keep the current sequence(s) open?")
        if keepfiles == None:
            return
        # Clear cache list if we are replacing files
        if keepfiles == False:
            clearSeqList()
        # If we're keeping files, cache old values from text box as a tuple unless sequence is empty
        if keepfiles == True:
            OldTxtSequence = txt_sequence.get(1.0, tk.END).strip()
            OldProtName = ent_protname.get().strip()
            if OldTxtSequence != "":
                sequence_cache_list[current_file_index] = (OldProtName,OldTxtSequence)
    # Prompt for Protein Name (numbers will be added)
    ProteinName = simpledialog.askstring(title='Protein Name',prompt='Protein Name:')
    if ProteinName == None: # i.e. if Cancel is pressed
        return
    with open(filepath, "r") as f:
        filecontents = f.readlines() # Gives list of lines
        # If .tsv file contains header line with titles, get rid of header
        firstline=filecontents[0]
        if 'Bracket' in firstline or 'Strategy' in firstline:
            filecontents.pop(0)
        # Loop through each line and find a text block that looks like a bracket
        bracket_re_match = r'[\s\"]*([\w\(\),>-]+\[.+\][\w\(\)\[\],>-]+)[\s\"]*' # Not pretty, but matches any SEGMENT[n]SEGMENT with or without annotations
        bracket_count=0
        for i,line in enumerate(filecontents):
            linesplit = re.split('\t',line)
            try:
                Sequence=re.search(bracket_re_match,line)[1]
                bracket_count+=1
            except:
                print(f'line {i+1} = {line}')
                report(f'File format error; no bracket found in line {i+1}')
                return
            # If we've gotten this far, we have a bracket; add it to our queue of sequences to load, as a tuple with protein name
            ProteinFullName=ProteinName+" - Line Number "+str(i+1)
            sequence_cache_list.append((ProteinFullName,Sequence))
            # First line only: load in text editor
            if i==0:
                total_file_number = len(sequence_cache_list)
                current_file_index = total_file_number-1
                var_currentfile.set(f'{current_file_index+1}')
                # Erase old values in text box and insert new ones
                txt_sequence.delete(1.0,tk.END)
                ent_protname.delete(0,tk.END)
                txt_sequence.insert(tk.END, Sequence)
                ent_protname.insert(tk.END, ProteinFullName)
                firstfilepath = filepath
            # Last file only: update global variables and report result
            if i==len(filecontents)-1:
                # Update total file number
                total_file_number = len(sequence_cache_list)
                # Update string variables in file navigation menu
                var_totalfiles.set(f'/{total_file_number}')
                # Load preview image of current open sequence (should be first in list)
                syntaxHighlight()
                loadPreviewImage()
                report(f'Loaded {bracket_count} brackets from {filepath}')

def btnOpenText(): # Open single .txt file using filedialogue
    filelist = filedialog.askopenfilenames(
        initialdir="C:/Users/MainFrame/Desktop/",
        title="Open Text File(s)",
        filetypes=(("Text Files (.txt, .fasta)", ".txt .fasta"),)
        )
    if not filelist:
        return
    global sequence_cache_list
    global current_file_index
    global total_file_number
    # Dialog popup; returns True to keep files, False to get rid of files
    # Ignore if only sequence is empty
    if len(sequence_cache_list)>0:
        keepfiles = messagebox.askyesnocancel("Question", "Do you want to keep the current sequence(s) open?")
        if keepfiles == None:
            return
        # Clear cache list if we are replacing files
        if keepfiles == False:
            clearSeqList()
        # If we're keeping files, cache old values from text box as a tuple unless sequence is empty
        if keepfiles == True:
            OldTxtSequence = txt_sequence.get(1.0, tk.END).strip()
            OldProtName = ent_protname.get().strip()
            if OldTxtSequence != "":
                sequence_cache_list[current_file_index] = (OldProtName,OldTxtSequence)
    for (i,filepath) in enumerate(filelist):
        with open(filepath, "r") as f:
            filecontents = f.readlines() # Gives list of lines
            # If .txt file contains header line beginning with >, interpret as protein name and the rest as sequence
            if filecontents[0][0] == '>':
                ProteinName = filecontents[0][1:].strip()
                Sequence = ''.join(filecontents[1:]).strip()
            # Otherwise interpret entire file contents as sequence, and keep protein name from before
            else:
                ProteinName = ""
                Sequence = ''.join(filecontents).strip()
            # Append onto our sequence cache list as a tuple
            sequence_cache_list.append((ProteinName,Sequence))
            # First file only: load in text editor
            if i==0:
                total_file_number = len(sequence_cache_list)
                current_file_index = total_file_number-1
                var_currentfile.set(f'{current_file_index+1}')
                # Erase old values in text box and insert new ones
                txt_sequence.delete(1.0,tk.END)
                ent_protname.delete(0,tk.END)
                txt_sequence.insert(tk.END, Sequence)
                ent_protname.insert(tk.END, ProteinName)
                firstfilepath = filepath
            # Last file only: update global variables and report result
            if i==len(filelist)-1:
                # Update total file number
                total_file_number = len(sequence_cache_list)
                # Update string variables in file navigation menu
                var_totalfiles.set(f'/{total_file_number}')
                # Load preview image of current open sequence (should be first in list)
                syntaxHighlight()
                loadPreviewImage()
                if len(filelist)==1:
                    report(f'Imported text file {filepath}')
                else:
                    report(f'Imported text file {firstfilepath} and {len(filelist)-1} others')

def btnSaveAll(): # Save all open sequences as individual text and/or image files (in specified formats)
    # If I make changes here, also change in btnSaveOne if relevant
    if reload_bracketmaker_on_run == True:
        importlib.reload(bracketmaker)
    global usr_save_txt
    global usr_save_svg
    global usr_save_png
    # First make sure text box isn't empty...BracketMaker can't run on an empty file
    ProteinName = ent_protname.get().strip()
    InputSequence = txt_sequence.get(1.0, tk.END).strip()
    if len(InputSequence) == 0:
        report("Error: Please enter Protein Sequence")
        return
    # Check format and exit if incorrect (only matters for images)
    if (highlightError() == False or checkImageSettings() == False) and (usr_save_png.get() == True or usr_save_svg.get() == True):
        return
    # Get list of file types
    filetypeslist = []
    if usr_save_txt.get() == True:
        filetypeslist.append(('Text File (.txt)','.txt'))
    if usr_save_svg.get() == True:
        filetypeslist.append(('Vector Graphics (.svg)','.svg'))
    if usr_save_png.get() == True:
        filetypeslist.append(('Raster Image (.png)','.png'))
    if len(filetypeslist) == 0:
        report('Error: Please select 1 or more file types to save as.')
        return
    # Ask folder to save in
    if len(filetypeslist) == 1:
        savetitle=f'Save folder for {filetypeslist[0][1]} files'
        extensionslist = f'({filetypeslist[0][1]})'
    else:
        extensionslist = (x[1] for x in filetypeslist)
        savetitle=f'Save folder for {extensionslist} files'
    # Prompt for save folder and file root name
    savefolder = filedialog.askdirectory(mustexist=True,title=savetitle)
    if savefolder == None:
        return
    filerootname = simpledialog.askstring(title=f'Save File Name {extensionslist}',prompt='Enter root file name (numbered extensions will automatically be added).')
    filename=f'{savefolder}/{filerootname}'
    # Cache current sequence to its position on the sequence list
    global sequence_cache_list
    global current_file_index
    OldTxtSequence = txt_sequence.get(1.0, tk.END).strip()
    OldProtName = ent_protname.get().strip()
    sequence_cache_list[current_file_index] = (OldProtName,OldTxtSequence)
    # Loop through every sequence currently open and run the appropriate functions to save as txt, png, and/or svg
    savedfilenames = []
    for i,seqtuple in enumerate(sequence_cache_list):
        report(f'Saving sequence {i+1} of {len(sequence_cache_list)}, please wait...')
        ProteinName = seqtuple[0]
        InputSequence = seqtuple[1]
        numberextension = str(i+1).zfill(len(str(len(sequence_cache_list)))) # Funky, but this auto detects how many zeroes to add
        fileheader = f'>{filerootname} strategy {numberextension}'
        # Save all file types in order
        if usr_save_txt.get() == True:
            txtfilename = f'{filename}_{numberextension}.txt'
            savedfilenames.append(txtfilename)
            with open(txtfilename, "w") as f:
                fileheader = ProteinName
                filecontents = InputSequence
                f.write(f'>{fileheader}\n{filecontents}')
        if usr_save_png.get() == True:
            pngfilename = f'{filename}_{numberextension}.png'
            savedfilenames.append(pngfilename)
            # Get preview image and save copy to this filename
            loadPreviewImage()
            shutil.copy('bm_cache/preview.png', pngfilename)
        if usr_save_svg.get() == True:
            svgfilename = f'{filename}_{numberextension}.svg'
            savedfilenames.append(svgfilename)
            # Save text box contents to a temporary file
            tempfilepath = 'bm_cache/currentsequence.txt'
            with open(tempfilepath, "w") as f:
                f.write(f'>{ProteinName}\n{InputSequence}')
            SVGFileContents = runMainProgram(InputSequence)
            with open(svgfilename, "w") as f:
                f.write(SVGFileContents)
    # Report names of saved files
    andothers = ''
    if len(savedfilenames) > 1:
        andothers = f'and {len(savedfilenames)-1} others'
    report(f'Saved file {savedfilenames[0]} {andothers}')

def btnSaveOne(): # Save file as specified formats.
    # Note - If changes are made to this fxn, also change in btnSaveAll if relevant
    if reload_bracketmaker_on_run == True:
        importlib.reload(bracketmaker)
    global usr_save_txt
    global usr_save_svg
    global usr_save_png
    # First make sure text box isn't empty...BracketMaker can't run on an empty file
    ProteinName = ent_protname.get().strip()
    InputSequence = txt_sequence.get(1.0, tk.END).strip()
    if len(InputSequence) == 0:
        report("Error: Please enter Protein Sequence")
        return
    # Check format and settings and exit if incorrect (only matters for images)
    if (highlightError() == False or checkImageSettings() == False) and (usr_save_png.get() == True or usr_save_svg.get() == True):
        return
    # Get list of file types
    filetypeslist = []
    if usr_save_txt.get() == True:
        filetypeslist.append(('Text File (.txt)','.txt'))
    if usr_save_svg.get() == True:
        filetypeslist.append(('Vector Graphics (.svg)','.svg'))
    if usr_save_png.get() == True:
        filetypeslist.append(('Raster Image (.png)','.png'))
    if len(filetypeslist) == 0:
        report('Error: Please select 1 or more file types to save as.')
        return
    if len(filetypeslist) == 1:
        savetitle='Save File'
    else:
        savetitle='Save Files'
    # Ask to save as file
    # This will provide a dropdown menu to select file types, but regardless of choice will use results from check boxes
    filepath = filedialog.asksaveasfile(
        mode='w',
        title = savetitle,
        defaultextension='.txt',
        filetypes=tuple(filetypeslist))
    # Get file name without any extension
    filename = os.path.splitext(filepath.name)[0]
    # Save all file types in order
    savedfilenames = []
    if usr_save_txt.get() == True:
        txtfilename = f'{filename}.txt'
        savedfilenames.append(txtfilename)
        with open(txtfilename, "w") as f:
            fileheader = ent_protname.get().strip()
            filecontents = txt_sequence.get(1.0, tk.END).strip()
            f.write(f'>{fileheader}\n{filecontents}')
    if usr_save_png.get() == True:
        pngfilename = f'{filename}.png'
        savedfilenames.append(pngfilename)
        # Get preview image and save copy to this filename
        loadPreviewImage()
        shutil.copy('bm_cache/preview.png', pngfilename)
    if usr_save_svg.get() == True:
        svgfilename = f'{filename}.svg'
        savedfilenames.append(svgfilename)
        # Save text box contents to a temporary file
        tempfilepath = 'bm_cache/currentsequence.txt'
        with open(tempfilepath, "w") as f:
            f.write(f'>{ProteinName}\n{InputSequence}')
        SVGFileContents = runMainProgram(InputSequence)
        with open(svgfilename, "w") as f:
            f.write(SVGFileContents)
        # Update global variable with last SVG (only used for quick re-run)
        global global_lastsvg
        global_lastsvg = svgfilename
    # Report names of saved files
    andothers = ''
    if len(savedfilenames) > 1:
        andothers = f'and {len(savedfilenames)-1} others'
    report(f'Saved file {filepath.name} {andothers}')
    return

def btnSaveSettings(): # Save current Image Settings to a custom .bmsettings file for later retrieval
    # Prompt to save a .bmsettings file
    filename = filedialog.asksaveasfilename(
        title = 'Save BracketMaker settings to file',
        defaultextension='.bmsettings',
        filetypes=(("Settings (.bmsettings)", ".bmsettings"),)
        )
    if not filename:
        return
    # Check format of image settings; cancel if incorrect
    if checkImageSettings()==False:
        return
    # Settings are saved as python code. When this file is read, it will execute the below statements and set variables.
    with open(filename,'w') as f:
        f.write(f'usr_align_segments_at_top.set({usr_align_segments_at_top.get()})\n')
        f.write(f'usr_show_thioester_placeholder.set({usr_show_thioester_placeholder.get()})\n')
        f.write(f'usr_thioester.set("{usr_thioester.get()}")\n')
        f.write(f'usr_show_segment_labels.set({usr_show_segment_labels.get()})\n')
        f.write(f'usr_show_special_aa_annotations.set({usr_show_special_aa_annotations.get()})\n')
        f.write(f'usr_colors = {tuple(usr_colors)}\n')
        f.write(f'usr_highlight_color = "{usr_highlight_color}"\n')
        f.write(f'usr_show_rxn_label_placeholder.set({usr_show_rxn_label_placeholder.get()})\n')
        f.write(f'usr_rxn_label_placeholder_text.set("{usr_rxn_label_placeholder_text.get()}")\n')
        f.write(f'usr_aa_start_number.set("{usr_aa_start_number.get()}")\n')
        f.write(f'usr_poor_thioesters.set("{usr_poor_thioesters.get()}")\n')
        f.write(f'usr_highlight_poor_thioesters.set({usr_highlight_poor_thioesters.get()})\n')
        f.write(f'usr_px_per_aa.set("{usr_px_per_aa.get()}")\n')
        f.write(f'usr_segment_spacing.set("{usr_segment_spacing.get()}")\n')
    report(f'Successfully saved settings file {filename}')
    return

def checkAutoBracketSettings():
    # Looks for values in Sort Settings boxes, returns True if correct, reports error and returns False if incorrect
    numberoflines=ent_number_of_brackets.get().strip()
    if not re.fullmatch('\d+', numberoflines):
        report(f'Please enter only numbers for Max Schemes')
        return False
    cyspg1 = ent_internal_cys_pg.get().strip()
    if cyspg1=='':
        report(f'Please enter Default Cys Protecting Group')
        return False
    if not re.fullmatch('[\w-]+',cyspg1):
        report(f'Please use only alphanumeric characters for Default Cys Protecting Group')
        return False
    cyspg2 = ent_nterm_cys_pg.get().strip()
    if cyspg2 == '':
        return True
    if not re.fullmatch('[\w-]+',cyspg2):
        report(f'Please use only alphanumeric characters for N-terminal Cys Protecting Group')
        return False
    return True

def checkImageSettings():
    # Validates values in settings boxes and reports an error if incorrect; returns True/False
    # AA start number must be a number
    try:
        x = int(usr_aa_start_number.get())
    except:
        report("Please enter only numbers for AA Start Number")
        return False
    # px values must be numbers
    try:
        x = int(usr_px_per_aa.get())
    except:
        report("Please enter only numbers for AA width (px)")
        return False
    try:
        x = int(usr_segment_spacing.get())
    except:
        report("Please enter only numbers for Segment Spacing (px)")
        return False
    return True

def clearPreviewImage(): # Remove everything from the image preview frame
    for widget in frm_imgpreview.winfo_children():
        widget.pack_forget()

def clearSeqList():
    # Removes all sequences from queue and starts over at 0
    global current_file_index
    global total_file_number
    global sequence_cache_list
    global var_currentfile
    global var_totalfiles
    sequence_cache_list = []
    current_file_index = 0
    total_file_number=1
    # Update string variables in file navigation menu
    var_currentfile.set(f'{current_file_index+1}')
    var_totalfiles.set(f'/{total_file_number}')


def displayBracketScores():
    # Change color of labels to grey
    lbl_maxpath.config(fg='grey')
    lbl_avgpath.config(fg='grey')
    lbl_steps.config(fg='grey')
    # Get input sequence and pass to scoring fxn with current user settings
    InputSequence = txt_sequence.get(1.0, tk.END).strip()
    if usr_penalize_poor_thioesters.get()==True:
        goodthioesterslist=getGoodThioesters(usr_poor_thioesters.get())
    else:
        goodthioesterslist=getGoodThioesters('')
    if usr_penalize_non_cys_thiols.get()==True:
        noncysthiolpenalty=2 # This is the default penalty in the scoring function
    else:
        noncysthiolpenalty=0
    # Check that format is correct, then get scores
    formatinfo = bracketmaker.checkFormat(InputSequence)
    if formatinfo[0]==True: # Must have at least one bracket in place
        ScoreDict=bracketmaker.scoreBracket(InputSequence,good_thioesters=goodthioesterslist,non_cys_thiol_penalty=noncysthiolpenalty)
        # Get the relevant scores and convert to strings, and set text variables...labels should automatically update
        var_maxpath.set(str(ScoreDict["maxpath"]))
        avgpath_float="{:.2f}".format(ScoreDict["avgpath"])
        var_avgpath.set(avgpath_float)
        var_steps.set(str(ScoreDict["steps"]))
        # Change color of labels to black
        lbl_maxpath.config(fg='black')
        lbl_avgpath.config(fg='black')
        lbl_steps.config(fg='black')

def getDarkerColor(InputTuple,fraction=0.5): # For automatically picking a dark Stroke color from a given Fill color, converts a hex to a darker color and returns hex code
    newr=int((InputTuple[0]*fraction))
    newg=int((InputTuple[1]*fraction))
    newb=int((InputTuple[2]*fraction))
    return "#{:02x}{:02x}{:02x}".format(newr,newg,newb)

def getGoodThioesters(InputString): # Converts a string (e.g. 'I,L,K,V,T') into a list of good thioesters by subtracting these from the list of all other possible characters. This is because all BracketMaker and Auto-Bracket functions use a list of "good" thioesters, but it is more intuitive for the user to specify poor ones.
    ThioestersList=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ.,-()[]')
    for Char in InputString:
        if Char in ThioestersList:
            ThioestersList.remove(Char)
    return ThioestersList

def getTextFileName(InputFilepath): # From filepath/folder/filename.txt, extract "filename"
    justfile=os.path.basename(InputFilepath)
    minusextension = os.path.splitext(justfile)[0]
    return minusextension

def highlightError(): # Detect syntax errors and highlight incorrect characters, displaying info at the bottom
    TxtSequenceContents = txt_sequence.get('1.0',tk.END).strip()
    formatinfo = bracketmaker.checkFormat(TxtSequenceContents) # Tuple containing info about format; see bracketmaker.py
    report('')
    if formatinfo[0]==False:
        if formatinfo[1] == 'invalid':
            texterror_context = formatinfo[2]
            texterror_index = formatinfo[3]
            texterror_char = formatinfo[4]
            report(f'Syntax error: encountered invalid character {texterror_char} within {texterror_context}.')
            # Delete last tag highlighting invalid character if found
            if txt_sequence.tag_names().count('invalidchar') != 0:
                txt_sequence.tag_delete('invalidchar')
            # Add tag highlighting invalid character
            txt_sequence.tag_add('invalidchar', f'1.{texterror_index}', f'1.{texterror_index+1}')
            txt_sequence.tag_config('invalidchar', background='yellow')
        elif formatinfo[1] == 'notags':
            report('Must have at least 1 ligation tag [#] to generate image (this will be fixed soon).')
        elif formatinfo[1] != 'list': # List of segments is the only valid format left
            report(f'Syntax error: report from bracketmaker = {formatinfo[1]}')
    return formatinfo[0]


def loadPreviewImage(force=False): # Load preview image - run the current sequence in BracketMaker using 'return_canvas = True', which returns Drawing object without text formatting. Then save as preview PNG, and display this PNG in the frame.
    # Refresh list of swatch buttons
    loadSwatchButtons()
    # Refresh bracket scores
    displayBracketScores()
    # Unless forced, this will only run if the global variable is active
    if (auto_preview_image.get() == True and show_preview_image.get()==True) or force==True:
        if reload_bracketmaker_on_run == True:
            importlib.reload(bracketmaker)
        # First make sure text box isn't empty...BracketMaker can't run on an empty file
        ProteinName = ent_protname.get().strip()
        InputSequenceBefore = txt_sequence.get(1.0, tk.END).strip()
        InputSequence = "".join(InputSequenceBefore.split()) # Remove ALL white space from sequence
        # To avoid errors in the BracketMaker error highlighting function, white space is not allowed in the sequence; replace the active text box with non-whitespace version
        if len(InputSequenceBefore)!=len(InputSequence):
            print(f'Before = {InputSequenceBefore} Length = {len(InputSequenceBefore)}\nAfter =  {InputSequence} Length = {len(InputSequence)}')
            txt_sequence.delete(1.0, tk.END)
            txt_sequence.insert(1.0, InputSequence)
            syntaxHighlight() # Re-do syntax highlighting
        # NOTE: Disabled requirement for protein name; it's fine if it's blank
        # if len(ProteinName) == 0:
        #     report("Error: Please enter Protein Name")
        #     return
        if len(InputSequence) == 0:
            report("Error: Please enter Protein Sequence")
            return
        # Check format and settings and exit if incorrect
        if highlightError() == False or checkImageSettings() == False:
            return # No need to report anything; other fxns will do that
        # Run BracketMaker in 'preview mode' using the text from text box:
        try:
            preview_canvas = runMainProgram(InputSequence,return_canvas=True)
        # This will fail if it is not in the correct format for BracketMaker, or if bracketmaker fails for any other reason
        except:
            print(f'BRACKETMAKER PREVIEW FAILURE')
            if force==True:
                report("Error: bracketmaker.makeBracketFigure failure")
            return
        # Cache image as PNG
        preview_canvas.savePng('bm_cache/preview.png')
        # Load the cached image
        render = resizedImage('bm_cache/preview.png')
        # Remove everything else from the preview frame
        clearPreviewImage()
        # Define label with preview image
        lbl_imgpreview.config(image=render)
        lbl_imgpreview.image = render # have to include this line for the image to show up
        lbl_imgpreview.pack()
        report('')

def loadPreviewImageFromKey(event): # Truncated version of above; calls this function on the given key press but does not send the key input through (initially mapped to Shift+Return)
    loadPreviewImage(force=True)
    return 'break'

def loadSequenceStrings(event,cacheindex=None,cache=True): # Function used by navigation menu and Open File buttons to load a new sequence
    # Check contents of entry box; if cache=True (i.e. when switching between open files), this will be preserved in the current location; be sure to pass cacheindex
    # If current text box contents are to be deleted/replaced, pass cache=False
    # This updates all string variables related to the total sequence cache, checks the global 'current position' and loads that sequence, generating preview image and syntax highlighting
    InputString = ent_filenumber.get()
    # Get global variables
    global total_file_number
    global current_file_index
    if cacheindex == None:
        cacheindex = current_file_index
    # Only allow numbers, otherwise set back
    if not re.fullmatch('\d+',InputString):
        report("Please enter a valid file number")
        var_currentfile.set(f'{current_file_index+1}')
        return
    # Set total file number variable
    var_totalfiles.set(f'/{total_file_number}')
    # Ensure number is valid; otherwise set to highest or lowest available
    InputNumber = int(InputString)
    # If too low, load sequence 1
    if InputNumber<1:
        var_currentfile.set('1')
        InputNumber = 1
    elif InputNumber>total_file_number:
        var_currentfile.set(f'{total_file_number}')
        InputNumber = total_file_number
    # If a valid number, cache the current text contents in place
    report("")
    global sequence_cache_list
    if cache==True:
        OldTxtSequence = txt_sequence.get(1.0, tk.END).strip()
        OldProtName = ent_protname.get().strip()
        sequence_cache_list[cacheindex] = (OldProtName,OldTxtSequence)
    # Load sequence cached at the requested index
    current_file_index = InputNumber-1
    ProteinName = sequence_cache_list[current_file_index][0]
    Sequence = sequence_cache_list[current_file_index][1]
    # Erase old values in text box and insert new ones
    txt_sequence.delete(1.0,tk.END)
    ent_protname.delete(0,tk.END)
    txt_sequence.insert(tk.END, Sequence)
    ent_protname.insert(tk.END, ProteinName)
    # Highlight syntax and generate preview image
    syntaxHighlight()
    loadPreviewImage()

def loadSwatchButtons():
    # Count number of segments in current text box by splitting at ligation tags
    SegmentsOnly=re.split(r'\[.*?\]', txt_sequence.get('1.0',tk.END))
    while '' in SegmentsOnly:
        SegmentsOnly.remove('') # Remove blank spaces in between explicit rxn steps
    shortcolorlist = usr_colors[0:len(SegmentsOnly)]
    # Get rid of previous buttons
    destroylist = []
    for child in frm_segment_colors.children.values():
        if isinstance(child,tk.Button):
            destroylist.append(child)
    for widget in destroylist:
        widget.destroy()
    # Loop through first x colors in current color list
    for i,pair in enumerate(shortcolorlist):
        fillhex = pair[0]
        strokehex = pair[1]
        # Create a button tied to this index, labeled with this color, within the frm_segment_colors frame
        tk.Button(frm_segment_colors,
            command=lambda x = i: btnChangeSegmentColor(x),
            image=pixel,
            width=20,
            height=20,
            relief='sunken',
            bg=fillhex).grid(row=0,column=i)
        # For Mac version, create a text label underneath these buttons
        tk.Label(frm_segment_colors,
            text=f'{i+1}',
            fg=fillhex,
            font=('Arial Bold',12)).grid(row=1,column=i)

def report(InputString): # Display text on the bottom of the screen
    str_footer.set(InputString)

def resizedImage(filename): # Open an image file and return an ImageTk object resized to match the max preview width
    load = Image.open(filename)
    # Resize if necessary
    orig_w=load.size[0]
    orig_h=load.size[1]
    if orig_w>img_preview_max_width or orig_h>img_preview_max_height:
        resize_w = img_preview_max_width/orig_w
        resize_h = img_preview_max_height/orig_h
        if resize_w<resize_h: # Pick which dimension to constrain the resize to
            resize_factor=resize_w
        else:
            resize_factor=resize_h
        load=load.resize((round(orig_w*resize_factor),round(orig_h*resize_factor)))
    # Render as an 'imageTk' object
    return ImageTk.PhotoImage(load)

def runMainProgram(InputSequence,return_canvas=False): # Execute bracketmaker.makeBracketFigure with current user settings in GUI
    # Set default thioester to blank if the setting is toggled False
    thioester_to_pass = usr_thioester.get()
    if usr_show_thioester_placeholder.get()==False:
        thioester_to_pass=''
    # Other settings are within the fxn call below
    return bracketmaker.makeBracketFigure(InputSequence,return_canvas=return_canvas,align_segments_at_top=usr_align_segments_at_top.get(),thioester=thioester_to_pass,show_segment_labels=usr_show_segment_labels.get(),show_special_aa_annotations=usr_show_special_aa_annotations.get(),colors=tuple(usr_colors),highlight_color=usr_highlight_color,show_rxn_label_placeholder=usr_show_rxn_label_placeholder.get(),rxn_label_placeholder_text=usr_rxn_label_placeholder_text.get(),aa_start_number=int(usr_aa_start_number.get()),good_thioesters=getGoodThioesters(usr_poor_thioesters.get()),highlight_poor_thioesters=usr_highlight_poor_thioesters.get(),highlight_non_cys_thiols=usr_highlight_non_cys_thiols.get(),px_per_aa=int(usr_px_per_aa.get()),segment_spacing=int(usr_segment_spacing.get()))

def syntaxHighlight(): # Highlight syntax with color-coding in sequence text box
    TxtSequenceContents = txt_sequence.get(1.0, tk.END).strip()
    # Remove all tags from current text
    for tag in txt_sequence.tag_names():
        txt_sequence.tag_delete(tag)
    # Initialize list and toggle switches
    TxtHighlightList = []
    ToggleNewTag = True # Triggers creation of a new tag upon encountering the next character
    InsideBracket = AfterHyphen = InsideParentheses = False # Bracket, hyphen, parentheses toggle
    ColorQueue = ['black'] # Only populated by default color at first
    InsideBracketColor = inside_bracket_color
    InsideParenthesesColor = usr_highlight_color
    AfterHyphenColor = after_hyphen_color
    i_max = len(TxtSequenceContents)
    isubtract = 0 # Number to subtract from index when starting new line
    LineNumber = 1
    # sub-fxn; add tag to list with a given start index using the ongoing color queue; it spans from the current index all the way to the end
    def addTag(startindex):
        if len(TxtHighlightList) > 0:
            TxtHighlightList[-1][1] = startindex # If there is a previous tag on the list, end it here
        TxtHighlightList.append([startindex,tk.END,ColorQueue[0]]) # Start new tag from here to end with the color in our queue, [start,end,color]
    # Start index counter
    for i in range(0,i_max): # Loop through positions
        Char = TxtSequenceContents[i]
        TagIndex = f'{LineNumber}.{i-isubtract}'
        # White space characters will change the index passed to new tags
        if re.fullmatch(r'\n',Char):
            LineNumber += 1
            isubtract = i+1 # Resets next to index 0 like a typewriter
        # Special characters will add colors to the queue and start creation of a new tag
        if Char == '[':
            InsideBracket = True # Toggle switch
            ColorQueue.insert(0,InsideBracketColor) # Add color to queue
            addTag(TagIndex) # End the last tag (if active) and start a new one
            # Opening bracket will close an 'after hyphen' tag as well
            if AfterHyphen == True:
                AfterHyphen = False
                ColorQueue.remove(AfterHyphenColor)
        elif Char == '(' and AfterHyphen == False:
            InsideParentheses = True # Toggle switch
            ColorQueue.insert(0,InsideParenthesesColor) # Add color to queue
            addTag(TagIndex) # End the last tag (if active) and start a new one
        elif Char == '-' and InsideParentheses == InsideBracket == AfterHyphen == False:
            AfterHyphen = True # Toggle switch
            ColorQueue.insert(0,AfterHyphenColor) # Add color to queue
            addTag(TagIndex) # End the last tag (if active) and start a new one
        # If we are 'inside' one of these tags and encounter a closing character, remove the color from the queue and trigger a new tag to start on next character
        elif InsideBracket == True and Char == ']':
            InsideBracket = False
            ColorQueue.remove(InsideBracketColor)
            ToggleNewTag = True
        elif InsideParentheses == True and Char == ')':
            InsideParentheses = False
            ColorQueue.remove(InsideParenthesesColor)
            ToggleNewTag = True
        # Finally, if this isn't a special character and the 'new tag' trigger is on, start a new tag with the last color in the queue
        elif ToggleNewTag == True:
            ToggleNewTag = False
            addTag(TagIndex)
    # Add tags
    for item in TxtHighlightList:
        tag_startindex = item[0]
        tag_endindex = item[1]
        tag_color = item[2]
        tag_name = str(item)
        txt_sequence.tag_add(tag_name, tag_startindex, tag_endindex) # Add tag with name,index1,index2
        txt_sequence.tag_config(tag_name, foreground=tag_color) # Define properties of tag

def syntaxHighlightAuto(event): # Automatically execute syntax highlighting and refresh preview image upon key press within text widget
    syntaxHighlight()
    loadPreviewImage()

def togglePenaltyThioester():
    # This fxn is triggered when changing check boxes in AutoBracket settings referring to poor-thioester and non-cys-thiol penalties
    # Set the Highlight X checkboxes to match these
    usr_highlight_poor_thioesters.set(usr_penalize_poor_thioesters.get())
    # This should automatically trigger a new loadPreviewImage due to tracing

def togglePenaltyThiol():
    usr_highlight_non_cys_thiols.set(usr_penalize_non_cys_thiols.get())

def toggleShowAutoBracket(event): # Expand or collapse the AutoBracket menu
    global show_autobracket
    # If turning on, add the frame to the canvas and update the string variable
    if show_autobracket == False:
        show_autobracket=True
        frm_autobracket.grid(row=4,columnspan=99,ipadx=5,ipady=5)
        str_autobracket.set('Auto-Bracket v')
        return
    # If turning off, remove the frame from the canvas and update the string variable
    show_autobracket=False
    str_autobracket.set('Auto-Bracket >')
    frm_autobracket.grid_forget()

def toggleShowImageSettings(event): # Expand or collapse the Figure Settings menu
    global show_image_settings
    # If turning on, add the frame to the canvas and update the string variable
    if show_image_settings == False:
        show_image_settings=True
        frm_image_settings.grid(row=4,columnspan=99,ipadx=5,ipady=5)
        str_settings_dropdown.set('Figure Settings v')
        return
    # If turning off, remove the frame from the canvas and update the string variable
    show_image_settings=False
    str_settings_dropdown.set('Figure Settings >')
    frm_image_settings.grid_forget()

def toggleShowPreview(): # Show or hide the preview image
    global show_preview_image
    global auto_preview_image
    if show_preview_image.get() == False: # If turning off...
        clearPreviewImage()
        # Turn off auto preview image and grey out text box
        auto_preview_image.set(False)
        s_auto_preview_image.config(state=tk.DISABLED)
    else: # If turning on instead
        s_auto_preview_image.config(state=tk.NORMAL)
        loadPreviewImage(force=True)


### --------------------------------
### III...GUI WINDOW ASSEMBLY
### --------------------------------

# "Traced" tkinter variables will call the specified fxn whenever they are changed
# Colors are not included in this list; they have their own function for loading values
usr_align_segments_at_top.trace('w',lambda *args: loadPreviewImage())
usr_show_thioester_placeholder.trace('w',lambda *args: loadPreviewImage())
usr_thioester.trace('w',lambda *args: loadPreviewImage())
usr_show_segment_labels.trace('w',lambda *args: loadPreviewImage())
usr_show_special_aa_annotations.trace('w',lambda *args: loadPreviewImage())
usr_show_rxn_label_placeholder.trace('w',lambda *args: loadPreviewImage())
usr_rxn_label_placeholder_text.trace('w',lambda *args: loadPreviewImage())
usr_aa_start_number.trace('w',lambda *args: loadPreviewImage())
usr_poor_thioesters.trace('w',lambda *args: loadPreviewImage())
usr_highlight_poor_thioesters.trace('w',lambda *args: loadPreviewImage())
usr_highlight_non_cys_thiols.trace('w',lambda *args: loadPreviewImage())
usr_px_per_aa.trace('w',lambda *args: loadPreviewImage())
usr_segment_spacing.trace('w',lambda *args: loadPreviewImage())
usr_penalize_poor_thioesters.trace('w',lambda *args: togglePenaltyThioester()) # Will load new preview image as well
usr_penalize_non_cys_thiols.trace('w',lambda *args: togglePenaltyThiol())
##
## MASTER FRAME - pack everything ultimately into this
##
# Frame encompassing entire GUI window
frm_master = tk.Frame()
frm_master.pack()


##
## HEADER TEXT - at top
##
frm_header = tk.Frame(frm_master)
frm_header.grid(row=0,columnspan=99)
# Welcome text at top of screen
lbl_welcome = tk.Label(text='Welcome to BracketMaker!\nEnter protein sequence below, or load a .txt or .fasta file.\nFor help building brackets, see the User Manual available on GitHub.',
    master=frm_header)
lbl_welcome.pack()


##
## MASTER FRAME FOR LEFT HALF OF SCREEN - TEXT EDITOR AND IMAGE PREVIEW
##
frm_master_left=tk.Frame(master=frm_master)
frm_master_left.grid(row=2,column=2,sticky='n')

# Frame encompassing Protein Name and Sequence boxes
frm_texteditor = tk.Frame(master=frm_master_left)
frm_texteditor.grid(row=2,column=2,padx=10,sticky='n')

# Protein name (single line)
frm_protname = tk.Frame(master=frm_texteditor)
frm_protname.grid(row=0,columnspan=99)
lbl_protname = tk.Label(master=frm_protname, text="Protein Name / File Header >")
lbl_protname.grid(row=0,column=0,sticky='e')
ent_protname = tk.Entry(master=frm_protname,width=90)
ent_protname.grid(row=0,column=1,sticky='w')

# Sequence box
txt_sequence = scrolledtext.ScrolledText(
    master=frm_texteditor,
    height=12)
# Bind syntax highlighting function
txt_sequence.bind('<KeyRelease>', syntaxHighlightAuto)
# Bind preview image function
txt_sequence.bind('<Shift-Return>',loadPreviewImageFromKey)
txt_sequence.bind('<space>',loadPreviewImageFromKey)
txt_sequence.bind('<Return>',loadPreviewImageFromKey)
txt_sequence.bind('<Tab>',loadPreviewImageFromKey)
txt_sequence.grid(row=2,columnspan=99)

# Protein name and sequence initial values
ent_protname.insert(tk.END, ent_protname_initial)
txt_sequence.insert(tk.END, txt_sequence_initial)
syntaxHighlight()

# NAVIGATION BUTTONS
# Frame within text editor encompassing file navigation buttons
frm_navigation = tk.Frame(master=frm_texteditor)
frm_navigation.grid(row=6,columnspan=99)

# Close button
btn_close = tk.Button(master=frm_navigation,
    text = 'Close (-)',
    fg='#882124',
    command=btnNavClose,
    width=10,
    height=1)
btn_close.grid(row=0,column=1)

# Prev arrow
btn_prevfile = tk.Button(master=frm_navigation,
    text="<--",
    command = btnNavPrevFile,
    width=18,
    height=1)
btn_prevfile.grid(row=0,column=2,padx=5)

# Current / Total sequences
ent_filenumber = tk.Entry(master=frm_navigation,
    width=3,
    textvariable=var_currentfile)
ent_filenumber.grid(row=0,column=4)
ent_filenumber.bind('<Return>',loadSequenceStrings)
ent_filenumber.bind('<FocusOut>',loadSequenceStrings)

lbl_totalfiles = tk.Label(master=frm_navigation,
    textvariable=var_totalfiles)
lbl_totalfiles.grid(row=0,column=6)

# Next arrow
btn_nextfile = tk.Button(master=frm_navigation,
    text="-->",
    command = btnNavNextFile,
    width=18,
    height=1)
btn_nextfile.grid(row=0,column=8,padx=5)

# Clone button
btn_clone = tk.Button(master=frm_navigation,
    text = 'Clone (+)',
    fg='#1f7d29',
    command=btnNavClone,
    width=10,
    height=1)
btn_clone.grid(row=0,column=10)

# IMAGE PREVIEW WINDOW
frm_imgpreview = tk.Frame(master=frm_master_left)
frm_imgpreview.grid(row=4,column=2)

lbl_imgpreview = tk.Label(master=frm_imgpreview)
if auto_load_last == True:
    render = resizedImage('bm_cache/preview.png')
    lbl_imgpreview = tk.Label(master=frm_imgpreview,image=render)
lbl_imgpreview.pack(expand=True,fill="both")


##
## MASTER FRAME FOR RIGHT SIDE - BUTTONS AND SETTINGS
##
frm_master_right=tk.Frame(master=frm_master)
frm_master_right.grid(row=2,column=4,sticky='n')

# OPEN AND SAVE BUTTONS, i.e. main buttons
# Frame encompassing main screen buttons
frm_mainbuttons = tk.Frame(master=frm_master_right)
frm_mainbuttons.grid(row=2)

btn_open = tk.Button(
    text="Open Text File(s)",
    command=btnOpenText,
    width=18,
    height=2,
    master=frm_mainbuttons)
btn_open.grid(row=0,column=0,padx=4,pady=4)

btn_openlist = tk.Button(
    text="Open From List",
    command=btnOpenList,
    width=18,
    height=2,
    master=frm_mainbuttons)
btn_openlist.grid(row=0,column=1,padx=4,pady=4)

# Save buttons with check boxes for filetypes
btn_saveone = tk.Button(
    text="Save Current As...",
    command=btnSaveOne,
    width=18,
    height=2,
    master=frm_mainbuttons)
btn_saveone.grid(row=1,column=0,padx=4,pady=4)

frm_savetypes = tk.Frame(frm_mainbuttons)
frm_savetypes.grid(row=1,rowspan=2,column=1)

s_save_txt = tk.Checkbutton(frm_savetypes,text='.txt (Plain Text)',variable=usr_save_txt)
s_save_txt.grid(row=0,sticky='w')
s_save_svg = tk.Checkbutton(frm_savetypes,text='.svg (Vector Graphics)',variable=usr_save_svg)
s_save_svg.grid(row=1,sticky='w')
s_save_png = tk.Checkbutton(frm_savetypes,text='.png (Preview PNG)',variable=usr_save_png)
s_save_png.grid(row=2,sticky='w')

##
## FIGURE SETTINGS WINDOW
##
frm_settings = tk.Frame(master=frm_master_right)
frm_settings.grid(row=3,sticky='n')
# Label with dropdown functionality
str_settings_dropdown=tk.StringVar(window)
str_settings_dropdown.set('Figure Settings >')
lbl_image_settings = tk.Label(master=frm_settings,textvariable=str_settings_dropdown,font=("Arial Bold",12),cursor='hand2')
lbl_image_settings.grid(row=3,columnspan=99)
lbl_image_settings.bind('<Button>',toggleShowImageSettings)
# Sub-frame containing all settings that will be passed to bracketmaker
frm_image_settings=tk.Frame(frm_settings,bd=3,relief=tk.RIDGE)

## Color swatch buttons
pixel = tk.PhotoImage(width=1,height=1) # virtual pixel with no info, useful for making square buttons
frm_segment_colors = tk.Frame(frm_image_settings)
frm_segment_colors.grid(row=1,columnspan=99)
# Populate list of color swatches; if text box is empty, this should only be one swatch
loadSwatchButtons()






# AA start number
lbl_aa_start_number = tk.Label(frm_image_settings, text='First AA Number:')
lbl_aa_start_number.grid(row=2,column=0,sticky='e')
ent_aa_start_number = tk.Entry(frm_image_settings,textvariable=usr_aa_start_number)
ent_aa_start_number.grid(row=2,column=1)

# Poor thioester and thiol highlights
chk_highlight_poor_thioesters = tk.Checkbutton(master=frm_image_settings,text='Highlight Poor Thioesters:', variable=usr_highlight_poor_thioesters,wraplength=150)
chk_highlight_poor_thioesters.grid(row=3,column=0,sticky='w')
ent_poor_thioesters = tk.Entry(master=frm_image_settings,textvariable=usr_poor_thioesters)
ent_poor_thioesters.grid(row=3,column=1,sticky='e')
chk_highlight_non_cys_thiols = tk.Checkbutton(master=frm_image_settings,text='Highlight Non-Cys Thiols', variable=usr_highlight_non_cys_thiols)
chk_highlight_non_cys_thiols.grid(row=4,column=0,columnspan=99)

# Toggle various elements of figure
lbl_toggle_text = tk.Label(frm_image_settings, text='Toggle text placeholders:', font=('Arial Bold',10))
lbl_toggle_text.grid(row=5,columnspan=99)
chk_show_segment_labels = tk.Checkbutton(frm_image_settings, text='Segment labels',variable=usr_show_segment_labels)
chk_show_segment_labels.grid(row=6,columnspan=99)
chk_show_special_aa_annotations = tk.Checkbutton(frm_image_settings, text='Special AA annotations',variable=usr_show_special_aa_annotations)
chk_show_special_aa_annotations.grid(row=7,columnspan=99)
chk_show_rxn_label_placeholder = tk.Checkbutton(frm_image_settings, text='Ligation reactions:',variable=usr_show_rxn_label_placeholder)
chk_show_rxn_label_placeholder.grid(row=8,column=0,sticky='e')
ent_rxn_label_placeholder_text = tk.Entry(frm_image_settings,textvariable=usr_rxn_label_placeholder_text)
ent_rxn_label_placeholder_text.grid(row=8,column=1,sticky='w')
chk_show_thioester_placeholder = tk.Checkbutton(frm_image_settings,text='C-term thioesters:',variable=usr_show_thioester_placeholder)
chk_show_thioester_placeholder.grid(row=9,column=0,sticky='e')
ent_thioester = tk.Entry(frm_image_settings,textvariable=usr_thioester)
ent_thioester.grid(row=9,column=1,sticky='w')

# Other image settings
lbl_other_image_settings = tk.Label(master=frm_image_settings,text='Other image settings:',font=('Arial Bold',10))
lbl_other_image_settings.grid(row=10,columnspan=99)
lbl_px_per_aa = tk.Label(master=frm_image_settings,text='AA width (px):')
lbl_px_per_aa.grid(row=11,column=0,sticky='e')
ent_px_per_aa = tk.Entry(master=frm_image_settings,textvariable=usr_px_per_aa)
ent_px_per_aa.grid(row=11,column=1,sticky='w')
lbl_segment_spacing = tk.Label(master=frm_image_settings,text='Segment spacing (px):')
lbl_segment_spacing.grid(row=12,column=0,sticky='e')
ent_segment_spacing = tk.Entry(master=frm_image_settings,textvariable=usr_segment_spacing)
ent_segment_spacing.grid(row=12,column=1,sticky='w')

# Save and Load Settings buttons
btn_save_settings = tk.Button(frm_image_settings,text='Save Settings',command=btnSaveSettings,height=1,width=12)
btn_save_settings.grid(row=15,column=0)
btn_load_settings = tk.Button(frm_image_settings,text='Load Settings',command=btnLoadSettings,height=1,width=12)
btn_load_settings.grid(row=15,column=1)

##
## AUTOBRACKET MENU
##
frm_autobracket_master=tk.Frame(frm_master_right) # This frame contains only the drop-down string and the frame with everything else inside
frm_autobracket_master.grid(row=4,sticky='n')
# Label with dropdown functionality
str_autobracket=tk.StringVar(window)
str_autobracket.set('Auto-Bracket >')
lbl_autobracket = tk.Label(master=frm_autobracket_master,textvariable=str_autobracket,font=("Arial Bold",12),cursor='hand2')
lbl_autobracket.grid(row=0,columnspan=99)
lbl_autobracket.bind('<Button>',toggleShowAutoBracket)

# Frame that is collapsed or expanded
frm_autobracket = tk.Frame(frm_autobracket_master,bd=3,relief=tk.RIDGE)

# Table of bracket scores
frm_scores = tk.Frame(frm_autobracket)
frm_scores.grid(row=1,column=0,columnspan=99)

lbl_score_header = tk.Label(frm_scores,text='Current Bracket Info:',font=("Arial Bold",10))
lbl_score_header.grid(row=1,column=0,columnspan=3)

lbl_maxpath_title = tk.Label(frm_scores,text='Max\nPath')
lbl_maxpath_title.grid(row=2,column=0)

lbl_avgpath_title = tk.Label(frm_scores,text='Avg\nPath')
lbl_avgpath_title.grid(row=2,column=1)

lbl_steps_title = tk.Label(frm_scores,text='Rxn\nSteps')
lbl_steps_title.grid(row=2,column=2)

lbl_maxpath = tk.Label(frm_scores,textvariable=var_maxpath)
lbl_maxpath.grid(row=3,column=0)

lbl_avgpath = tk.Label(frm_scores,textvariable=var_avgpath)
lbl_avgpath.grid(row=3,column=1)

lbl_steps = tk.Label(frm_scores,textvariable=var_steps)
lbl_steps.grid(row=3,column=2)

# AutoBracket functionality
lbl_autobracket_help = tk.Label(frm_autobracket,text='To auto-fill, insert [] between segments and click\nBuild Bracket, or import an Aligator file.')
lbl_autobracket_help.grid(row=10,column=0,columnspan=99)

btn_autobracket = tk.Button(
    text="Build Bracket(s)",
    command=btnAutoBracket,
    width=18,
    height=2,
    master=frm_autobracket)
btn_autobracket.grid(row=11,column=0,padx=2,pady=2)

btn_aligator = tk.Button(
    text="Aligator Import",
    command=btnAligatorImport,
    width=18,
    height=2,
    master=frm_autobracket)
btn_aligator.grid(row=11,column=1,padx=2,pady=2)

chk_save_autobracket_file = tk.Checkbutton(frm_autobracket,text='Save AutoBracket Results (.tsv)',variable=usr_save_autobracket_file)
chk_save_autobracket_file.grid(row=13,column=0,columnspan=99,sticky='n')

lbl_number_of_brackets = tk.Label(master=frm_autobracket,
    text='How many different possible brackets? (0 for all)',wraplength=150)
lbl_number_of_brackets.grid(row=14,column=0,sticky='e')
ent_number_of_brackets = tk.Entry(master=frm_autobracket)
ent_number_of_brackets.insert(tk.END,'1')
ent_number_of_brackets.grid(row=14,column=1,sticky='w')

lbl_autobracket_settings = tk.Label(frm_autobracket,text='AutoBracket Settings',font=("Arial Bold",10))
lbl_autobracket_settings.grid(row=15,column=0,columnspan=99,sticky='n')

# Sub-frame including autobracket settings (with their defaults)
frm_autobracket_settings = tk.Frame(master=frm_autobracket)
frm_autobracket_settings.grid(row=16,columnspan=3,padx=4,pady=4)

chk_penalize_thioesters = tk.Checkbutton(frm_autobracket_settings,text='Penalize Poor Thioesters (+1 Path)',variable=usr_penalize_poor_thioesters)
chk_penalize_thioesters.grid(row=0,columnspan=99,sticky='n')

chk_penalize_non_cys_thiols = tk.Checkbutton(frm_autobracket_settings,text='Penalize Non-Cys Thiols (+2 Path)',variable=usr_penalize_non_cys_thiols)
chk_penalize_non_cys_thiols.grid(row=1,columnspan=99,sticky='n')

lbl_internal_cys_pg = tk.Label(master=frm_autobracket_settings,
    text='Cys Protecting Group:',wraplength=150)
lbl_internal_cys_pg.grid(row=3,column=0,sticky='e')
ent_internal_cys_pg = tk.Entry(master=frm_autobracket_settings)
ent_internal_cys_pg.insert(tk.END,'Acm')
ent_internal_cys_pg.grid(row=3,column=1,sticky='w')

lbl_nterm_cys_pg = tk.Label(master=frm_autobracket_settings,
    text='N-terminal Cys Protecting Group (if different):',wraplength=150)
lbl_nterm_cys_pg.grid(row=4,column=0,sticky='e')
ent_nterm_cys_pg = tk.Entry(master=frm_autobracket_settings)
ent_nterm_cys_pg.insert(tk.END,'Tfa-Thz')
ent_nterm_cys_pg.grid(row=4,column=1,sticky='w')

# Sub-frame for one-pot toggle check boxes
frm_one_pot_toggles = tk.Frame(frm_autobracket)
frm_one_pot_toggles.grid(row=17,columnspan=3,padx=4,pady=4)

lbl_one_pot_steps = tk.Label(master=frm_one_pot_toggles,text='Toggle one-pot steps:',wraplength=50)
lbl_one_pot_steps.grid(rowspan=99,column=0,sticky='e')

chk_one_pot_desulfurization = tk.Checkbutton(frm_one_pot_toggles,text='Desulfurization',variable=usr_one_pot_desulfurization)
chk_one_pot_desulfurization.grid(row=0,column=1,sticky='w')
chk_one_pot_internal_cys = tk.Checkbutton(frm_one_pot_toggles,text='Internal Cys Deprotection',variable=usr_one_pot_internal_cys)
chk_one_pot_internal_cys.grid(row=1,column=1,sticky='w')
chk_one_pot_nterm_cys = tk.Checkbutton(frm_one_pot_toggles,text='N-Terminal Cys Deprotection',variable=usr_one_potnterm_cys)
chk_one_pot_nterm_cys.grid(row=2,column=1,sticky='w')

##
## FOOTER CONTENT
##
frm_footer = tk.Frame(master=frm_master)
frm_footer.grid(row=99,column=0,columnspan=99)

str_footer = tk.StringVar() # To make changes, update this string variable
str_footer.set("")

lbl_footer = tk.Label(master=frm_footer,
    textvariable=str_footer,
    height=1,
    anchor="w",
    justify=tk.LEFT)
lbl_footer.grid(sticky='w')


###
### IV...RUN GUI
###
# Print welcome text upon successful launch
print("Welcome to BracketMaker!\nCreated by Judah Evangelista, Kay Lab, University of Utah\nNow running BracketMaker GUI...")

# Run the main loop of the GUI
if open_gui == True:
    window.mainloop()
