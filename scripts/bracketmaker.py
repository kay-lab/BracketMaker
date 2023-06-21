#! /usr/bin/env python

#-------------Welcome to BracketMaker (work in progress)------------
#------------------Judah Evangelista, Kay Lab, University of Utah---
#-------------------------------------------------------------------
#-------------------------------------------------------------------


#-------------------------------------------------------------------
#--------------HEADER - Imported Modules and Dev Stuff--------------
#-------------------------------------------------------------------
import re
import sys
import os
from PIL import ImageColor # Has some functions for converting between hex/RGB color values
import drawSvg as draw
from collections import Counter
import itertools
import datetime # For naming files


# For dev - report time to screen
import time
start_time = time.time()
def gettime(reporttext):
    if __name__ == '__main__':
        print(f'{reporttext} -  {time.time()-start_time} seconds')
###########

#-------------------------------------------------------------------
#--------------DEFAULT VARIABLES ACCESSED BY FUNCTIONS BELOW--------
#-------------------------------------------------------------------

default_good_thioesters=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ.,()[]') # Keep ) for special characters; assume they're good


#-------------------------------------------------------------------
#--------------FUNCTION DEFS - alphabetical order-------------------
#-------------------------------------------------------------------

def autoBracket(InputSequence,returnnumber=1000,internal_cys_pg='Acm',n_term_cys_pg='Tfa-Thz',good_thioesters=default_good_thioesters,non_cys_thiol_penalty=2,one_pot_cys_deprotection=False,one_pot_desulfurization=False,one_pot_nterm_cys_deprotection=True,return_scores=False,output_scores_file=False,output_partial_brackets=False,filename='autobracketresults',quickmode=False,verbose=False,desulf_pairs={"A":"C","V":"Pen"},record_loop_times=False):
    fxn_call_time=time.time()
    # Takes input in the form of a partially-finished bracket, with blank brackets [] in desired ligation sites
    # Makes an output .tsv with all strategy scores if output_scores_file is True; regardless, returns a list of final strategies sorted in score order
    # INTIAL PROCESSING - break into segment and tag list
    if verbose==True:
        print(f'Running autoBracket on {InputSequence}')
    SegmentLigationList = re.split(r'(\[.*?\])', InputSequence)
    # If there are any blank spaces, this is a space between two explicit rxn steps; reformat to interpret as single tag
    TempList=[]
    BlankSpaceFound=False
    for item in SegmentLigationList:
        if BlankSpaceFound==False and item!='':
            TempList.append(item)
            LastItem=item
        if BlankSpaceFound==True:
            TempList[-1]=LastItem+item
            BlankSpaceFound=False
        if item=='':
            BlankSpaceFound=True
    SegmentLigationList=TempList[:]
    # Break into list of segments and ligation tags
    SegmentList = SegmentLigationList[0::2]
    TagList = SegmentLigationList[1::2]
    # Determine which thioester positions are GOOD and which are BAD; used for scoring/sorting tags before placement
    global GoodLigationsList
    GoodLigationsList = []
    for Segment in SegmentList[0:len(SegmentList)-1]:
        RightAA = Segment[-1]
        if (RightAA in good_thioesters):
            GoodLigationsList.append(True)
        else:
            GoodLigationsList.append(False)
    # Initialize global strategy queue, final strategy list, and dictionaries with scores
    global StrategyQueue
    StrategyQueue = [tuple(SegmentLigationList)] # Enclosed within double brackets; a list of lists, currently with only one sub-list
    global FinalBracketList
    FinalBracketList=[]
    global BracketScoreDict
    BracketScoreDict={}
    # Begin processing strategies in the queue
    # This loop will end when either the final strategy list is sufficiently full, or there are no more strategies to process
    loops_processed=0
    process_start_time=time.time()-fxn_call_time
    LoopTimeDict={0:process_start_time} # Record the time to finish each loop
    DoneChecking=False
    NextLoopIsFinal=False
    while len(StrategyQueue)>0 and DoneChecking==False:
        if verbose==True:
            gettime(f'-----\nPLACING ALL [{loops_processed+1}] TAGS')
        processStrategyQueue(internal_cys_pg=internal_cys_pg,n_term_cys_pg=n_term_cys_pg,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty,one_pot_cys_deprotection=one_pot_cys_deprotection,one_pot_desulfurization=one_pot_desulfurization,one_pot_nterm_cys_deprotection=one_pot_nterm_cys_deprotection,desulf_pairs=desulf_pairs)
        loops_processed+=1
        # After each loop, if there are finished strategies, record the time at which this loop finished
        if len(FinalBracketList)>0:
            LoopTimeDict[loops_processed]=time.time()-fxn_call_time
        # QUICK MODE - if toggled, function will run 1 more loop after any strategies have been found, then end the search
        # In >99.99% of Aligator strategies from E. coli proteome test set, the best possible bracket is among the 1st or 2nd 'batch' of complete brackets
        if quickmode==True:
            if NextLoopIsFinal==True:
                DoneChecking=True
            if len(FinalBracketList)>=returnnumber:
                NextLoopIsFinal=True
        # If toggled to output partial brackets, the processed Strategy Queue (i.e., all possible placements of the next level number) will be reported as a file, similar to the final output list
        if output_partial_brackets==True:
            if verbose==True:
                gettime(f'-----\nDONE WITH LOOP {loops_processed}...writing partial brackets to file')
            with open(f'output/{filename}_{loops_processed}.tsv','w') as f:
                f.write('Strategy\tMin Yield\tAvg Yield\tMax Path\tSum Path\tAvg Path\tSteps\tThioester Penalty\tNumber of Tags\n')
                for BracketListForm in StrategyQueue:
                    Bracket=''.join(BracketListForm)
                    score = scoreBracket(Bracket,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty) # Dictionary containing all scores
                    f.write(f'{Bracket}\t{stringPercent(score["minyield"])}\t{stringPercent(score["avgyield"])}\t{score["maxpath"]}\t{score["sumpath"]}\t{score["avgpath"]}\t{score["steps"]}\t{score["thioester"]}\t{countTags(BracketListForm)}\n')
    # SORT FINAL STRATEGY LIST
    if verbose==True:
        gettime(f'Found {len(FinalBracketList)} possible strategies')
    SortedFinalBracketList = sorted(FinalBracketList, key = lambda x: (BracketScoreDict[x]["maxpath"], BracketScoreDict[x]["sumpath"],BracketScoreDict[x]["steps"])) # Sorts from lowest up, which is what we want with path length; make negative if using yields
    # Trim list to top N options
    if returnnumber<len(SortedFinalBracketList):
        FinalBracketList = SortedFinalBracketList[0:returnnumber]
    else:
        FinalBracketList = SortedFinalBracketList[:]
    if verbose==True:
        gettime(f'Sorted list to top {len(FinalBracketList)}')
    # Optional - create TSV with all reported scores
    if output_scores_file==True:
        with open(f'output/{filename}.tsv','w') as f:
            f.write('Strategy\tMin Yield\tAvg Yield\tMax Path\tSum Path\tAvg Path\tSteps\tThioester Penalty\n')
            for Bracket in FinalBracketList:
                score = BracketScoreDict[Bracket] # Dictionary containing all scores
                f.write(f'{Bracket}\t{stringPercent(score["minyield"])}\t{stringPercent(score["avgyield"])}\t{score["maxpath"]}\t{score["sumpath"]}\t{score["avgpath"]}\t{score["steps"]}\t{score["thioester"]}\n')
    # If toggled, return list of formatted and sorted brackets along with score dictionary
    if return_scores==True: # True when this fxn is called by Aligator import; False when called on its own
        BracketScoreList=[]
        for Bracket in FinalBracketList:
            BracketScoreTuple=(Bracket,BracketScoreDict[Bracket])
            BracketScoreList.append(BracketScoreTuple)
        if record_loop_times==True:
            return BracketScoreList,LoopTimeDict
        else:
            return BracketScoreList
    # Otherwise, return list of formatted and sorted brackets only
    return FinalBracketList


def checkFormat(InputString):
# Analyzes format of string and determines if it is the correct format for BracketMaker.
# If good, returns a tuple starting with True and the sequence information to pass to other functions
# If not good, returns a tuple with False, the type of error, and other info relevant to the error type
    InputString = InputString.strip()
    # First determine if this is a single string or multiple
    StringList = InputString.split()
    if len(StringList) > 1: # Multiple segments separated by white space
        # If this is a list of strings, check if any of them have ligation tags
        HasTags = False
        for x in StringList:
            if re.search(r'\[.*?\]', x):
                HasTags = True
        # If any segment has a ligation tag, ignore all white space and treat this as a single string
        if HasTags == True:
            InputString = ''.join(StringList)
        else:
            return (False,'list',StringList) # Otherwise, return list of segments
    # If this is a single string, detect any ligation tags []
    SegmentLigationList = re.split(r'(\[.*?\])', InputString)
    # At this point, fxn will only continue if this appears to be a string with at least one ligation tag
    # FIRST PASS CHECK FOR VALID CHARACTERS
    InsideParentheses = InsideBracket = AfterHyphen = False
    remember_bracket_i = remember_paren_i = 0
    charset_segment = r'[a-zA-Z\[\(-]' # Any letter, opening bracket, opening parentheses, and hyphen
    charset_parentheses = r'[\w\-\>\)]' # Any letter or number, hyphen -, right arrow >, closing parentheses
    charset_bracket = r'[\w\-\>\],\*]' # Any letter or number, hyphen -, right arrow >, closing bracket, comma, asterisk*
    charset_hyphen = r'[\w\(\)\-\[]' # Any letter or number, parentheses, hyphen, opening bracket
    # Loop through all characters in the string and check if they are valid given the context
    for i in range(0,len(InputString)):
        Char = InputString[i]
        # Check if it matches the current context, return an informative error if not
        if InsideBracket == True:
            if re.fullmatch(charset_bracket, Char) == None:
                return (False,'invalid','bracket',i,Char) # False, error type, context, character index, character
        elif InsideParentheses == True:
            if re.fullmatch(charset_parentheses, Char) == None:
                return (False,'invalid','parentheses',i,Char) # False, error type, context, character index, character
        elif AfterHyphen == True:
            if re.fullmatch(charset_hyphen, Char) == None:
                return (False,'invalid','C-terminus',i,Char) # False, error type, context, character index, character
        else:
            if re.fullmatch(charset_segment, Char) == None:
                return (False,'invalid','segment',i,Char) # False, error type, context, character index, character
        # Certain characters will trigger a toggle in character set for the next loop
        if Char == '-' and InsideBracket==InsideParentheses==False:
            AfterHyphen = True
        if Char == '[':
            InsideBracket = True
            AfterHyphen = False
            remember_bracket_char = Char
            remember_bracket_i = i
        if Char == '(':
            InsideParentheses = True
            remember_paren_char = Char
            remember_paren_i = i
        if Char == ']':
            InsideBracket = False
        if Char == ')':
            InsideParentheses = False
    # One last check - if it still thinks we're inside parentheses or a bracket, the opening bracket must be an invalid character
    if InsideBracket==True:
        return(False,'invalid','open bracket',remember_bracket_i,'[')
    if InsideParentheses==True:
        return(False,'invalid','open parentheses',remember_paren_i,'(')
    # If it passed all the checks, return a True tuple with variables to pass to the next functions
    return (True,InputString,SegmentLigationList)

def condenseAnnotation(InputList):
    # For labels above segments denoting special internal AA; combines two text labels into one
    # Input is a list of tuples of the form [(AAnumber, Label, CenterPoint),...]
    # Output is a single tuple of the same format...
    # Combine text into a list which will be condensed back into a string
    CombinedStringList=[]
    Labels=[x[1] for x in InputList]
    LastText=""
    for Label in Labels:
        # Label should contain text and a number; split into these parts
        SplitLabel=re.search('(.*?)(\d+)', Label)
        Text=SplitLabel[1]
        Number=SplitLabel[2]
        # If this text is the same as the last text, add a comma and just the number
        if Text==LastText:
            CombinedStringList.append(",")
            CombinedStringList.append(Number)
        # If this text is different from the last text, add a semicolon (except for first in list) unless numbers are subscripted, then the whole tag text and number
        if Text!=LastText:
            if len(CombinedStringList)>0 and preview_mode==True:
                CombinedStringList.append(";")
            CombinedStringList.append(Text)
            CombinedStringList.append(Number)
            # Finally, if this text was different from the last text, update it as the last text for the next loop
            LastText=Text
    CombinedString=''.join(CombinedStringList)
    # Get average center position
    CenterPositions = [x[2] for x in InputList]
    AvgPosition = sum(CenterPositions)/len(CenterPositions)
    # AAnumber is not meaningful, so just use the first AA number on the first item on the list
    FirstNumber=InputList[0][0]
    return (FirstNumber,CombinedString,AvgPosition)


def countTags(InputTuple):
    # Counts the number of tags that have been placed in a partially-built strategy
    LigationTagSearchPattern = r"\[.+\]" #Use + instead of *; this ignores empty tags
    t = 0
    for x in InputTuple:
        if re.fullmatch(LigationTagSearchPattern, x):
            t+=1
    return t


def drawBracket(InputVar, Level, level_shift=0, parent_id=0, lr="", good_thioesters=default_good_thioesters, changed_aas_list=[]):
    # This recursive fxn takes a Segment input in the form of a range (referenced to global sequence list) or a recursive range (multiple segments)
    # First, it draws a box for the Segment, placing it at the correct location depending on the segment's level, length, and position within the sequence.
    # Then, if this segment has internal ligations (i.e., is a list), it draws a bracket-style rxn line above leading to the L and R halves.
    # Finally, if this segment has internal ligations, run the fxn again on the L and R halves.
    desulf_pairs={"A":"C","V":"Pen"} # Define AAs that are changed upon "Desulfurization" keyword, After:Before
    # Get segment ID from global ID counter, then update counter
    global SegmentIDCounter
    SegmentID = SegmentIDCounter
    SegmentIDCounter += 1
    # Create new entry in segment dictionary with Segment ID as key; value is a sub-dictionary that will be filled with parameters for placing objects on the canvas; for now, it only contains the parent ID
    SegmentDict.update({
        SegmentID:{'parent':parent_id}
    })
    # Update parent segment dictionary entry to include L and R half IDs if this is a L or R half segment; also update this dictionary entry to include parent
    if lr == "l":
        SegmentDict[parent_id].update({"leftid":SegmentID})
    if lr == "r":
        SegmentDict[parent_id].update({"rightid":SegmentID})
    if lr == "c":
        SegmentDict[parent_id].update({"centerid":SegmentID})
    # TERMINAL SEGMENTS ONLY - Set this condition, which will be used later
    isTerminalSegment=False
    if isinstance(InputVar, range):
        isTerminalSegment = True
    # GET SEGMENT BOUNDARIES
    # Get boundaries of L and R segments
    FirstAAIndex = getLeftBound(InputVar)
    LastAAIndex = getRightBound(InputVar)
    SegmentLength = LastAAIndex-FirstAAIndex+1
    # GET SEGMENT LABEL FROM SEQUENCE LIST
    # Initial variables
    AnnotationDict={} # key = AA number, value = text of label (protecting group only for N/C-term AA, full AA identity for middle ones)
    HighlightFirstAA=False
    HighlightLastAA=False
    # GET FIRST AND LAST AA and their current identity in this segment; label version will be abbreviated if special AA
    FirstAA = SequenceList[FirstAAIndex]
    FirstAANumber = FirstAAIndex+aa_start_number_global
    LastAA = SequenceList[LastAAIndex]
    LastAANumber = LastAAIndex+aa_start_number_global
    # Update these in our dictionary
    SegmentDict[SegmentID].update({
        "firstAA":FirstAA,
        "firstAAnumber":FirstAANumber,
        "lastAA":LastAA,
        "lastAAnumber":LastAANumber
        })
    # Check if first AA is a special AA, and if so determine its current identity at this level
    if isinstance(FirstAA,dict):
        CheckNumber=Level
        if (FirstAAIndex,Level) in changed_aas_list:
            CheckNumber=Level+1
        CurrentFirstAA=FirstAA[CheckNumber]
        # If there is a protecting group (indicated by hyphen), add this to our annotations
        HyphenSplit=re.split("-", CurrentFirstAA, maxsplit=1)
        FirstAAShort=HyphenSplit[0]
        SegmentDict[SegmentID].update({"firstAA":FirstAAShort}) # unprotected AA identity will show on label
        if len(HyphenSplit)>1:
            FirstPG=HyphenSplit[1]
            AnnotationDict.update({FirstAANumber:FirstPG}) # Annotation text only includes PG
        # For desulfurizations and other AA changes, highlight the AA if it is currently not its final identity
        FinalFirstAA=FirstAA[BottomLevelNumber]
        FinalFirstAAShort=re.split("-",FinalFirstAA, maxsplit=1)[0]
        if FirstAAShort!=FinalFirstAAShort:
            HighlightFirstAA=True
    SegmentDict[SegmentID].update({"highlightfirstAA":HighlightFirstAA})
    # Do the same thing for last AA
    if isinstance(LastAA,dict):
        CheckNumber=Level
        if (LastAAIndex,Level) in changed_aas_list:
            CheckNumber=Level+1
        CurrentLastAA=LastAA[CheckNumber]
        # If there is a protecting group (indicated by hyphen), add this to our annotations
        HyphenSplit=re.split("-", CurrentLastAA, maxsplit=1)
        LastAAShort=HyphenSplit[0]
        SegmentDict[SegmentID].update({"lastAA":LastAAShort}) # unprotected AA identity will show on label
        if len(HyphenSplit)>1:
            LastPG=HyphenSplit[1]
            AnnotationDict.update({LastAANumber:LastPG}) # Annotation text only includes PG
        # For desulfurizations and other AA changes, highlight the AA if it is currently not its final identity
        FinalLastAA=LastAA[BottomLevelNumber]
        FinalLastAAShort=re.split("-",FinalLastAA, maxsplit=1)[0]
        if LastAAShort!=FinalLastAAShort:
            HighlightLastAA=True
    SegmentDict[SegmentID].update({"highlightlastAA":HighlightLastAA})
    # SPECIAL MIDDLE AAS...Show if different from final sequence, or on level of reaction
    for i in range(FirstAAIndex+1,LastAAIndex):
        ThisAA=SequenceList[i]
        if isinstance(ThisAA,dict):
            CheckNumber = Level
            # For single-segment rxns, if this is the level of rxn, treat the label as if it is one level higher
            if (i,Level) in changed_aas_list:
                CheckNumber = Level+1
            CurrentMiddleAA=ThisAA[CheckNumber]
            FinalMiddleAA=ThisAA[BottomLevelNumber]
            # If this is different from the final AA, or is equal to the final AA but the level above is not, annotate this amino acid above segment
            AnnotateThisAA=False
            if CurrentMiddleAA!=FinalMiddleAA:
                AnnotateThisAA=True
            if CheckNumber+1 in ThisAA:
                HigherMiddleAA=ThisAA[CheckNumber+1]
                if CurrentMiddleAA!=HigherMiddleAA and CurrentMiddleAA==FinalMiddleAA:
                    AnnotateThisAA=True
            if AnnotateThisAA==True:
                AnnotationDict.update({i+aa_start_number_global:formatPG(CurrentMiddleAA)})
    # Save our annotations dictionary; this will be used to build segment annotation tags later
    SegmentDict[SegmentID].update({"annotations":AnnotationDict})

    # GET SEGMENT DIMENSIONS AND CANVAS POSITION
    # Vertical spacing between segments is determined by rxn line length and segment height
    segment_y_spacing = line_vertical_length + segment_height/2
    # Default dimensions are based on sequence length, and are all relative to the bottom segment
    bottom_segment_x_left = -bottom_segment_width/2
    # X SHIFT - spacing left or right to give a 'funnel' shape to diagram, calculated based on SAME OR LOWER levels to left or right
    segment_x_left = bottom_segment_x_left + relativeCanvasPosition(FirstAAIndex) + getXShift(FirstAAIndex,LastAAIndex,Level)
    segment_y_bottom = (Level-BottomLevelNumber+level_shift)*(segment_height+segment_y_spacing)
    # segment_width = relativeCanvasPosition(SegmentLength)
    segment_width = px_per_aa_global*SegmentLength
    try:
        TopLevelNumber=max(AllLevelNumbers)+1 # 'true' top level of starting segments is 1 above top rxn
    except:
        TopLevelNumber=1 # For single segments
    DistanceFromTopLevel=TopLevelNumber-Level
    if DistanceFromTopLevel>0:
        segment_width=segment_width*(segment_funnel_scaling**DistanceFromTopLevel)
    # Other useful x,y coordinates (in px)
    segment_x_right = segment_x_left+segment_width
    segment_x_center = segment_x_left + segment_width/2
    segment_y_center = segment_y_bottom + segment_height/2
    segment_y_top = segment_y_bottom + segment_height
    # GET OTHER SEGMENT SHAPE PARAMETERS
    # COLORS - get relevant colors for segment based on AA bounds
    FillList = []
    StrokeList = []
    FoundNTermColor=False # Capture the first colors identified
    for i,fill,stroke in ColorStopList:
        if FirstAAIndex<=i<LastAAIndex:
            GradientPosition = (i-FirstAAIndex)/SegmentLength
            FillList.append((GradientPosition,fill))
            StrokeList.append((GradientPosition,stroke))
            # Update N-term fill color the first time a color match is found
            if FoundNTermColor==False:
                NTermFillColor=stroke
                NTermStrokeColor=stroke
                FoundNTermColor=True
            # Update C-term fill color every time a new color match is found
            CTermFillColor = stroke
            CTermStrokeColor = stroke
    # Pass some colors for N- and C-term text labels to our dictionary
    SegmentDict[SegmentID].update({
        "ctermfill":CTermFillColor,
        "ctermstroke":CTermStrokeColor,
        "ntermfill":NTermFillColor,
        "ntermstroke":NTermStrokeColor
        })
    # If single item in list, pass along a single hex string instead of a gradient list
    if len(FillList)==1:
        FillList = FillList[0][1]
        StrokeList = StrokeList[0][1]
    # DIFFERENCES FOR TERMINAL SEGMENTS
    fill_opacity = 1.0
    SegmentOpacity=fill_opacity
    AddStrokeWidth=0
    if isTerminalSegment==True:
        AddStrokeWidth=2
        SegmentOpacity=1.0
    # Define shape in Segment Dictionary entry
    SegmentDict[SegmentID].update({
        "x":segment_x_left,
        "y":segment_y_bottom,
        "width":segment_width,
        "center":segment_x_center,
        "height":segment_height,
        "fill":FillList,
        "strokewidth":segment_stroke_width+AddStrokeWidth,
        "stroke":StrokeList,
        "opacity":SegmentOpacity,
        "terminal":isTerminalSegment,
        "labelx":segment_x_center,
        "labely":segment_y_center-(segment_label_font_size/4)
        })
    # C TERM LABEL
    # This block of code draws the chemical group specified on the C terminus:
    CTermLabel=CTermLabelDict[LastAANumber]
    CTermWidth=getCharWidth(CTermLabel,annotation_font_size,sub_numbers=not preview_mode)
    if preview_mode==False:
        CTermLabel=subNumbers(CTermLabel)
    SegmentDict[SegmentID].update({
        "cterm":CTermLabel,
        "ctermx":segment_x_right+segment_stroke_width,
        "ctermwidth":CTermWidth
    })
    #---end DRAW CURRENT SEGMENT--#
    # DRAW REACTION LINE
    # MIDDLE SEGMENTS - Detect whether to draw a single line for an explicit rxn step, or a bracket line for a NCL step
    if isTerminalSegment==False:
        LigationTag = InputVar[1]
        StepsForLater = []
        explicit_step_found=False
        KeywordList = []
        for term in reversed(LigationTag):
            try:
                x = int(Keyword)
            except:
                if explicit_step_found==True:
                    StepsForLater.append(term)
                else:
                    KeywordList.append(term.rstrip('*'))
                    if '*' in term:
                        explicit_step_found=True
        KeywordList.reverse()
        StepsForLater.reverse()
        # Desulfurization keyword; if keyword contains any variant of desulf*, change it to 'Desulfurization' for proper interpretation
        TempList=[]
        for Keyword in KeywordList:
            if re.match('DESULF', Keyword.upper()):
                TempList.append('Desulfurization')
            else:
                TempList.append(Keyword)
        KeywordList=TempList[:]
        # If this segment undergoes a separate single-segment rxn, draw a single reaction line instead of bracket, and call function again with a similar inputvar, updated height.
        if explicit_step_found==True:
            # Only one X coordinate and two Y coordinates needed to draw rxn line
            rxnline_x_center = segment_x_left + segment_width/2 # Horizontal center of bottom segment
            rxnline_y_bottom = segment_y_bottom + segment_height # Top of bottom segment
            rxnline_y_top = rxnline_y_bottom + segment_y_spacing
            rxnline_y_center = (rxnline_y_bottom+rxnline_y_top)/2 # Used for placing text
            # GET REACTION LABEL PARAMETERS
            # Get rxn label text from keyword list
            RxnLabelText=KeywordList[0]
            if len(KeywordList)>1:
                for Keyword in KeywordList[1:]:
                    RxnLabelText+=f', {Keyword}'
            RxnLabelTextFull = f'{RxnLabelText}'
            RxnLabelX = rxnline_x_center + rxn_label_padding
            RxnLabelY = rxnline_y_center
            # Add line coordinates and label to segment dictionary entry
            SegmentDict[SegmentID].update({
                "singlelinex":rxnline_x_center,
                "singlelineyt":rxnline_y_top,
                "singlelineyb":rxnline_y_bottom,
                "rxnlabel":RxnLabelTextFull,
                "rxnlabelx":RxnLabelX,
                "rxnlabely":RxnLabelY
            })
            # Determine if any AA in bounds are affected by the current keyword, and will need to be interpreted differently in the next function:
            ChangedAAs = changed_aas_list # AAs changed at this level will be added to the list
            # Desulfurization keyword - add some extra keywords to the 'list to check' based on desulfurization pairs dictionary
            KeywordListToCheck=KeywordList[:]
            if 'Desulfurization' in KeywordList:
                for after_key in desulf_pairs:
                    before_key=desulf_pairs[after_key]
                    SubString=f'{before_key}>{after_key}'
                    KeywordListToCheck.append(SubString)
            # Check any AAs in these segments that are affected by a keyword change at some point
            for aai,aa in enumerate(SequenceList[FirstAAIndex:LastAAIndex+1]):
                if isinstance(aa,dict):
                    # AA changes such as C>A; build string matching keyword
                    checkstring1=checkstring2=None
                    if aa[Level+1]!=aa[Level]:
                        checkstring1 = aa[Level+1]+'>'+aa[Level]
                    tempsplit = aa[Level+1].split('-')
                    if len(tempsplit)>1 and tempsplit[0]==aa[Level]:
                        checkstring2 = '-'.join(tempsplit[1:])
                    if (checkstring1 in KeywordListToCheck) or (checkstring2 in KeywordListToCheck):
                        ChangedAAs.append((FirstAAIndex+aai,Level))
            # Input variable is a specific 3-item list. See makeNestedList() for details.
            # Construct a different version of the original InputVar, identical to this one except with the explicit step keyword removed
            LeftHalf = InputVar[0]
            RightHalf = InputVar[2]
            NewNestedList = [LeftHalf, StepsForLater, RightHalf]
            ls = level_shift+1
            drawBracket(NewNestedList,Level,parent_id=SegmentID,level_shift=ls,changed_aas_list=ChangedAAs,lr='c', good_thioesters=good_thioesters)
        # LIGATION JUNCTION - DRAW RXN LINE AND RUN FUNCTION ON BOTH HALVES
        else:
            # Input variable is a specific 3-item list. See makeNestedList() for details.
            LeftHalf = InputVar[0]
            LigationType = InputVar[1]
            RightHalf = InputVar[2]
            # Get AA bounds of left and right half segments
            l1 = FirstAAIndex
            l2 = getRightBound(LeftHalf)
            r1 = getLeftBound(RightHalf)
            r2 = LastAAIndex
            # Determine if this is a good or bad thioester; this will optionally highlight poor thioesters with different colored line
            ThioesterAA = SequenceList[l2]
            if isinstance(ThioesterAA,dict):
                ThioesterAA = ThioesterAA[Level]
            IsGoodThioester=False
            if ThioesterAA in good_thioesters:
                IsGoodThioester=True
            # Determine if this is a Cys or Ala thiol; others will be optionally highlighted
            ThiolAA=SequenceList[r1]
            if isinstance(ThiolAA,dict):
                ThiolAA = ThiolAA[Level]
            IsCysThiol=False
            if ThiolAA=='C' or ThiolAA=='A' or re.fullmatch('C-.*',ThiolAA):
                IsCysThiol=True
            # SIX PAIRS OF (x,y) COORDINATES ARE NEEDED TO DRAW THIS BRACKET
            # But only 3 unique values each for x and y
            rxnline_x_center = segment_x_left + segment_width/2 # Horizontal center of bottom segment
            rxnline_x_left = bottom_segment_x_left + relativeCanvasPosition((l1+l2)/2) + getXShift(l1,l2,Level+1) # Horizontal center of LeftHalf
            rxnline_x_right = bottom_segment_x_left + relativeCanvasPosition((r1+r2)/2) + getXShift(r1,r2,Level+1) # Horizontal center of RightHalf
            rxnline_y_bottom = segment_y_bottom + segment_height # Top of bottom segment
            rxnline_y_center = rxnline_y_bottom + line_vertical_length # Defined in user input section RXN LINE STYLE
            rxnline_y_top = rxnline_y_bottom + segment_y_spacing
            # Adjust rxn line center so it remains between left and right boundaries
            if rxnline_x_center < rxnline_x_left+rxn_label_padding:
                rxnline_x_center = rxnline_x_left+rxn_label_padding
            elif rxnline_x_center > rxnline_x_right-rxn_label_padding:
                rxnline_x_center = rxnline_x_right-rxn_label_padding
            # GET REACTION LABEL PARAMETERS
            # This adds a placeholder rxn label, if active
            RxnLabelText = rxn_label_placeholder_global
            RxnLabelX = rxnline_x_center + rxn_label_padding
            RxnLabelY = (rxnline_y_center+rxnline_y_bottom)/2
            # Add other reactions such as protecting group removal, desulfurization, etc. to placeholder text
            for keyword in KeywordList:
                try:
                    x = int(keyword)
                except:
                    RxnLabelText+=f', {keyword}'
            # Construct full string for rxn label text
            RxnLabelTextFull = f'{RxnLabelText}'
            # Add line coordinates and label to segment dictionary entry
            SegmentDict[SegmentID].update({
                "linexl":rxnline_x_left,
                "linexc":rxnline_x_center,
                "linexr":rxnline_x_right,
                "lineyt":rxnline_y_top,
                "lineyc":rxnline_y_center,
                "lineyb":rxnline_y_bottom,
                "rxnlabel":RxnLabelTextFull,
                "rxnlabelx":RxnLabelX,
                "rxnlabely":RxnLabelY,
                "goodthioester":IsGoodThioester,
                "cysthiol":IsCysThiol
            })
            #---begin RECURSION---#
            # IF INTERNAL LIGATIONS PRESENT, RUN THE FXN ON LEFT HALF, THEN RIGHT HALF
            drawBracket(LeftHalf,Level+1,parent_id=SegmentID,lr="l",level_shift=level_shift,changed_aas_list=changed_aas_list, good_thioesters=good_thioesters)
            drawBracket(RightHalf,Level+1,parent_id=SegmentID,lr="r",level_shift=level_shift,changed_aas_list=changed_aas_list, good_thioesters=good_thioesters)
            #---end RECURSION---#

def formatPG(InputString):
    # Converts PG format 'C-Acm' to 'C(Acm)'
    x = InputString.split('-')
    if len(x)==1:
        return x[0]
    else:
        return f'{x[0]}({"-".join(x[1:])})'

def getCharWidth(InputString,FontSize=12,sub_numbers=False):
    # This function determines the width of a text string in Arial font using measured widths for each individual character, and outputs the size in pixels (assuming 1:1 pixel ratio)
    # Used to determine minimum segment width based on the size of the label
    # Define width of all characters (measured from document using 120-point Arial bold font)
    widthdict = {
        33.4: split(r'.,'),
        34: split(r'ijlI'),
        40: split(r'()tf-'),
        47: split(r'r'),
        60: split(r'z1'),
        67: split(r'aceksvxyJ023456789'),
        74: split(r'bdghnopquFLTZ'),
        81: split(r'EPSVXY'),
        87: split(r'ABCDHKNRU'),
        94: split(r'wGOQ'),
        100: split(r'M'),
        107: split(r'm'),
        113: split(r'W')
    }
    # Default width for characters not in dictionary
    default_width = 70
    # Initialize loop
    TotalTextWidth = 0
    # Loop through all characters in the input string
    for Char in InputString:
        CharWidth = default_width
        # Compare to dictionary of width/character values and update to correct width
        for n in widthdict:
            if Char in widthdict[n]:
                CharWidth = n
        # Above values are based on 120-point font; adjust for actual font size. If number subscripting is on they will be smaller.
        if Char in ['0','1','2','3','4','5','6','7','8','9',','] and sub_numbers==True:
            ActualCharWidth = CharWidth/120*FontSize*10/12
        else:
            ActualCharWidth = CharWidth/120*FontSize
        TotalTextWidth += ActualCharWidth
    return TotalTextWidth

def getLeftBound(InputRange): # Returns the first number in a range or a recursive range
    if isinstance(InputRange,range):
        return InputRange[0]
    else:
        return getLeftBound(InputRange[0])

def getMaxPathLength(InputSequence,good_thioesters=default_good_thioesters,non_cys_thiol_penalty=2):
    # Runs below fxn but only returns max path length
    AllPathLengths = getPathLengths(InputSequence,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty).values()
    MaxPathLength = max(AllPathLengths)
    return MaxPathLength

def getPathLengths(InputSequence,PrevPathLength=0,Level=1,good_thioesters=default_good_thioesters,non_cys_thiol_penalty=2):
    # Recursive sub-fxn of scoreBracket
    # Input is a formatted BracketMaker strategy. Output is a dictionary of all the segments and their path lengths.
    global SegmentPathLengthDict
    # Define list of OK thiols that are exempt from non-Cys penalty. For now this is only C and A.
    ok_thiols=['C','A']
    # Split into list of segments & ligations
    LigationTagSearchPattern = r"(\[.*?\])"
    SegmentLigationList = re.split(LigationTagSearchPattern,InputSequence)
    # Initialize segment dictionary first time through fxn
    if Level==1:
        SegmentPathLengthDict = {}
    # Initialize other variables
    CurrentPathLength = PrevPathLength
    # If this is a single segment, no ligations, we have reached the end of a path. Update dictionary with max path length and end fxn
    if len(SegmentLigationList)==1:
        Segment = SegmentLigationList[0]
        # For scoring during autoBracket- Empty spaces where a bracket has not been placed are marked with a special character $
        if '$' in Segment:
            SegmentList=Segment.split('$')
        else:
            SegmentList=[Segment]
        for s in SegmentList:
            # Duplicate segments - if there is already a matching segment in dictionary, add a special character to denote this one
            while s in SegmentPathLengthDict:
                s=s+'+'
            SegmentPathLengthDict.update({s:CurrentPathLength})
        return
    # Loop through tag positions and identify the tag of the current level (there should only be one if formatted correctly!)
    Found=False
    for i in range(1,len(SegmentLigationList),2):
        if Found==False:
            LigationTag = SegmentLigationList[i]
            # Split into list of terms
            bracketandcomma = r'[\[\],]'
            LigationTagTerms = re.split(bracketandcomma,LigationTag)
            # Remove white space at beginning and end of this match
            while '' in LigationTagTerms:
                LigationTagTerms.remove('')
            # Find the tag matching the current level, and determine the total number of explicit (*) rxns in this tag, adding to the total path length so far
            if int(LigationTagTerms[0]) == Level:
                SplitPosition=i
                for Keyword in LigationTagTerms:
                    if Keyword.isnumeric()==True or Keyword[-1]=='*':
                        CurrentPathLength+=1
                Found=True
    # Split into left and right halves
    LeftList = SegmentLigationList[0:SplitPosition]
    RightList = SegmentLigationList[SplitPosition+1:]
    LeftSequence=''.join(LeftList)
    RightSequence=''.join(RightList)
    # Determine if this is a poor thioester, in which case path length is increased by one
    ThioesterAA = LeftSequence[-1]
    # If this is a special AA, use the identity of original AA
    if ThioesterAA==")":
        SplitByParen=re.split(r'[\(\)]', LeftSequence)
        while '' in SplitByParen:
            SplitByParen.remove('')
        # Use > and - to determine original AA
        if '>' in ThioesterAA:
            ThioesterAA=re.split('>',ThioesterAA)[0]
        if '-' in ThioesterAA:
            ThioesterAA=re.split('-',ThioesterAA)[0]
    if not ThioesterAA in good_thioesters:
        CurrentPathLength+=1
    # Determine if this is a poor thiol, in which case path length is increased by specified penalty
    # First, determine identity of thiol AA; if this is a special AA, get its identity before transformation
    ThiolAA=RightSequence[0]
    if ThiolAA=="(":
        FirstAAString=re.split(r'[\(\)]',RightSequence)[1]
        if '>' in FirstAAString:
            FirstAAString=re.split('>',FirstAAString)[0]
        if '-' in FirstAAString:
            FirstAAString=re.split('-',FirstAAString)[0]
        ThiolAA=FirstAAString
    # Penalize thiol AA if it is in dictionary
    if not ThiolAA in ok_thiols:
        CurrentPathLength+=non_cys_thiol_penalty
    # Run function again on left and right halves
    NextLevel = Level+1
    # RECURSION
    getPathLengths(LeftSequence,PrevPathLength=CurrentPathLength,Level=NextLevel,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty)
    getPathLengths(RightSequence,PrevPathLength=CurrentPathLength,Level=NextLevel,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty)
    # END RECURSION
    return SegmentPathLengthDict

def getRightBound(InputRange): # Returns the last number in a range or a recursive range
    if isinstance(InputRange,range):
        return InputRange[-1]
    else:
        return getRightBound(InputRange[2])

def getTagLevel(InputTagList,CheckLevel=1):
    # Gives the lowest tag number that can be placed in a valid empty site within a given list of tags
    LeftList=[]
    RightList=[]
    FoundTag=False
    for tag in InputTagList:
        if FoundTag==True:
            RightList.append(tag)
        else:
            if re.match(f'\[{CheckLevel}',tag):
                FoundTag=True
            else:
                LeftList.append(tag)
    if FoundTag==False:
        return CheckLevel
    if FoundTag==True:
        return min(getTagLevel(LeftList,CheckLevel+1),getTagLevel(RightList,CheckLevel+1))


# This function is called within drawBracket() to get the horizontal spacing of a given segment - needed to place segment box and reaction lines
def getXShift(ileft,iright,level): # ileft and iright are the index numbers of first and last amino acids
    for irange,(a,b) in enumerate(SequenceRangeList):
        if ileft==a:
            isegleft = irange
        if iright+1==b:
            isegright = irange
    LeftLevels = AllLevelNumbers[0:isegleft]
    if isegright==len(AllLevelNumbers):
        RightLevels=[]
    else:
        RightLevels=AllLevelNumbers[isegright:] # Avoids error for segments containing rightmost segment
    Count=0 # Positive if counted on left (shift segment right); negative if counted on right (shift segment left)
    for l in LeftLevels:
        if l<=level:
            Count+=1
    for r in RightLevels:
        if r<=level:
            Count-=1
    return Count*segment_horizontal_spacing


def getAligatorBrackets(InputFileContents,brackets_per_set=1,lines_to_read=10000,output_file=False,output_filename="AligatorBrackets",sort_by_bracket_score=True,internal_cys_pg='Acm',n_term_cys_pg='Tfa-Thz',good_thioesters=default_good_thioesters,one_pot_cys_deprotection=False,one_pot_desulfurization=False,one_pot_nterm_cys_deprotection=False,verbose=False,use_bracket_families=False,quickmode=False,non_cys_thiol_penalty=2):
    # Import Aligator function, used to generate the #1 bracket for each Aligator output strategy (i.e., list of segments) starting from the contents of an "All Strategies.txt" file
    if verbose==True:
        print(f'getAligatorBrackets was successfully called! Now reading file....')
    # If first line is header, pop it from the list
    if re.match('Strat', InputFileContents[0]):
        InputFileContents.pop(0)
    # Trim file contents if longer than specified read number
    if len(InputFileContents)>lines_to_read:
        InputFileContents=InputFileContents[0:lines_to_read]
    # Loop through Aligator lines and make brackets for strategies
    AligatorBracketList=[]
    AligatorLineNumber=1
    SingleSegmentList=[]
    BracketFamilyDict={} # For use_bracket_families mode
    BracketToKeyDict={}
    for lineindex,AligatorLine in enumerate(InputFileContents):
        if verbose==True:
            print(f'Reading line {lineindex} = {AligatorLine}')
        LineSplit=re.split('\s+', AligatorLine)
        # Store score for later
        AligatorScore=LineSplit.pop(0)
        # Remove blank spaces from list (incorrectly processed white-space characters)
        while '' in LineSplit:
            LineSplit.remove('')
        # Count number of segments
        NumberOfSegments=len(LineSplit)
        # Single segments will throw an error in autoBracket; remove these and report n/a later
        if NumberOfSegments==1:
            SingleSegmentList.append((LineSplit[0],AligatorLineNumber,AligatorScore))
        else:
            # BRACKET FAMILY MODE; Recommended for processing multiple large proteins, such as the top-1000 Aligator strategies from a large (>10 segment) protein.
            # If toggled, the program keeps track of every set of ligation junctions in a dictionary. If a similar segment set is encountered later in the list (same type & linear order of ligation junctions), the output list of best brackets will be the same as already calculated. The first of each 'family' is kept and duplicates are not processed, but are added to the family count.
            BracketFamilyKey=[]
            if use_bracket_families==True:
                # Begin building BracketFamilyKey, a list of tuples representing all ligation junctions in this set
                # Will be converted to a tuple of tuples when stored in bracket family dict
                for ileft in range(0,len(LineSplit)-1):
                    iright=ileft+1
                    LeftSegment=LineSplit[ileft]
                    RightSegment=LineSplit[iright]
                    # Tuple for each ligation junction is (IsGoodThioester, HasInternalCys, "ThiolAA") and will look something like (True, False, "C")
                    # Determine if this is a good thioester
                    ThioesterAA=LeftSegment[-1]
                    IsGoodThioester=False
                    if ThioesterAA in good_thioesters:
                        IsGoodThioester=True
                    # Determine if this ligation junction has any Cys in its segments (for determining if there may need to be a protection step for desulfurization, can be a difference between bracket families in rare cases)
                    HasInternalCys=False
                    if "C" in LeftSegment or "C" in RightSegment[1:]:
                        HasInternalCys=True
                    # Determine the thiol AA
                    ThiolAA=RightSegment[0]
                    # Add this tuple to list of ligation junctions
                    BracketFamilyKey.append((IsGoodThioester,HasInternalCys,ThiolAA))
                BracketFamilyKey=tuple(BracketFamilyKey)
            # First, try looking this up in the bracket family dict and increasing the count; if bracket family mode is off, this will automatically fail
            try:
                BracketFamilyDict[BracketFamilyKey]+=1
            # If bracket family mode is off, or if this family is new, run autoBracket on this line
            except:
                if verbose==True:
                    print(f'Calculating best bracket for line {AligatorLineNumber}, {NumberOfSegments} segments...')
                # LineSplit from above is just a list of segments; convert to a blank bracket
                StartBracket='[]'.join(LineSplit)
                # Run AutoBracket, and return the best bracket(s) for this line
                AutoBracketResults=autoBracket(StartBracket, internal_cys_pg=internal_cys_pg, n_term_cys_pg=n_term_cys_pg, good_thioesters=good_thioesters, output_scores_file=False, one_pot_cys_deprotection=one_pot_cys_deprotection, one_pot_desulfurization=one_pot_desulfurization, one_pot_nterm_cys_deprotection=one_pot_nterm_cys_deprotection, return_scores=True, returnnumber=brackets_per_set, record_loop_times=True, quickmode=quickmode,non_cys_thiol_penalty=non_cys_thiol_penalty)
                BracketsAndScores=AutoBracketResults[0]
                LoopTimes=AutoBracketResults[1]
                # Trim output list to only top bracket(s)
                if len(BracketsAndScores)>brackets_per_set:
                    BracketsAndScores=BracketsAndScores[0:brackets_per_set]
                for i,(Bracket,Scores) in enumerate(BracketsAndScores):
                    # Add a few more scores
                    Scores['aligator']=AligatorScore
                    Scores['a']=AligatorLineNumber
                    Scores['b']=i
                    Scores['segments']=NumberOfSegments
                    # Add scores for loop times
                    RecordedLoops=sorted(list(LoopTimes.keys()))
                    Scores['firstdoneloop']=RecordedLoops[1]
                    Scores['lastloop']=RecordedLoops[-1]
                    Scores['totaltime']=LoopTimes[RecordedLoops[-1]]
                    # Get ID and time of winning loop, by determining highest tag in bracket
                    AllLevelNumbers=re.findall(r'\d+',Bracket)
                    HighestLevelFound=int(max(AllLevelNumbers))
                    Scores['winningloop']=HighestLevelFound
                    Scores['winninglooptime']=LoopTimes[HighestLevelFound]
                    # Add to our master list of output strategies
                    AligatorBracketList.append((Bracket,Scores))
                # Bracket family mode - add an entry for this in our family dictionary, and (simple version) simply keep track of how many are in this family
                if use_bracket_families==True:
                    BracketFamilyDict.update({BracketFamilyKey:1})
                    # To look up this score later, save a different dictionary where the key is this calculated best bracket, and value is the bracket family key
                    BracketToKeyDict.update({BracketsAndScores[0][0]:BracketFamilyKey})
            # End loop through this Aligator line; increase line number for next loop
            AligatorLineNumber+=1
    # Done looping through Aligator lines
    # Print results to screen
    print(f'{len(BracketFamilyDict)} calculated brackets')
    # Sort list by scores of output brackets, grouped by number of segments, then by max path, then by sumpath, then by number of steps
    if sort_by_bracket_score==True:
        AligatorBracketList.sort(key=lambda x: (x[1]["segments"],x[1]["maxpath"],x[1]["sumpath"],x[1]["steps"])) # Normal order will have lower values first, which is what we want
    # OUTPUT TO FILE
    if output_file==True:
        with open(f'{output_filename}.tsv','w') as f:
            # Write file header
            f.write(f'AligatorRank\tBracketRank\tSHORTBRACKET\tNumber Segments\tAligator Score\tMax Path\tSum Path\tSteps\tAvg Path\tFULLBRACKET\tFamily Size\tFirst Finished Loop\tWinning Loop\tTotal Loops\tWinning Loop Time\tTotal Done Time\n')
            for BracketScoreTuple in AligatorBracketList:
                Bracket=BracketScoreTuple[0]
                scores=BracketScoreTuple[1]
                # Bracket family mode - report the number of similar segment sets encountered
                FamilyCount=1
                if Bracket in BracketToKeyDict:
                    BracketFamilyKey=BracketToKeyDict[Bracket]
                    FamilyCount=BracketFamilyDict[BracketFamilyKey]
                # Write bracket and scores to file
                f.write(f'{scores["a"]}\t{scores["b"]+1}\t{shortBracket(Bracket)}\t{scores["segments"]}\t{scores["aligator"]}\t{scores["maxpath"]}\t{scores["sumpath"]}\t{scores["steps"]}\t{round(scores["sumpath"]/scores["segments"],2)}\t{Bracket}\t{FamilyCount}\t{scores["firstdoneloop"]}\t{scores["winningloop"]}\t{scores["lastloop"]}\t{scores["winninglooptime"]}\t{scores["totaltime"]}\n')
            for Segment,AligatorLineNumber,AligatorScore in SingleSegmentList:
                f.write(f'{AligatorLineNumber}\tn/a\tn/a\t1\t{AligatorScore}\tn/a\tn/a\tn/a\tn/a\t{Segment}\n')
    # Return list of brackets and their scores
    return AligatorBracketList


def shortBracket(InputString):
    # Converts a full bracket PEPTIDE[2]PEPTIDE[1]PEPTIDE into a condensed tag format, [2][1]
    TagList = re.split(r'(\[.*?\])',InputString)
    ShortBracket=''.join(TagList[1::2])
    return ShortBracket


# Determines which level a keyword is removed by searching the ligation tag list
def keywordSearch(keyword,taglist,i): # i is the index of the segment
    # Determine what is the highest tag to start searching. This will be the highest tag immediately adjacent to the segment.
    leftnumber = 0
    if i>0:
        try:
            leftnumber = int(stringToList(taglist[i-1])[0])
        except:
            leftnumber = 0
    rightnumber = 0
    if i<len(taglist):
        try:
            rightnumber = int(stringToList(taglist[i])[0])
        except:
            rightnumber = 0
    top = max(leftnumber,rightnumber)
    for level in range(top,0,-1):
        rightboundfound = leftboundfound = False
        righttagtocheck = '[0]'
        lefttagtocheck = '[0]'
        if i<len(taglist):
            righttaglist = taglist[i:]
            for tag in righttaglist:
                try:
                    tagnumber = int(stringToList(tag)[0])
                except:
                    tagnumber=top+100
                if rightboundfound==False and tagnumber<=level:
                    rightboundfound=True
                    righttagtocheck = tag
        if i>0:
            lefttaglist = reversed(taglist[0:i])
            for tag in lefttaglist:
                try:
                    tagnumber = int(stringToList(tag)[0])
                except:
                    tagnumber=top+100
                if leftboundfound==False and tagnumber<=level:
                    leftboundfound=True
                    lefttagtocheck = tag
        for tag in righttagtocheck,lefttagtocheck:
            if int(stringToList(tag)[0])==level:
                for term in stringToList(tag):
                    if re.match(keyword,term):
                        return level
                    # Desulfurization keyword; perform a softer search for related terms
                    if re.match('DESULF',keyword.upper()) and re.match('DESULF',term.upper()):
                        return level
    # If fxn has not returned by now, there are no matches
    return 0


def listToString(InputList):
    # Converts a list of strings into a string that looks like a list, minus quote marks (e.g."[NCL,0]")
    OngoingString='['
    for i,s in enumerate(InputList):
        OngoingString+=str(s)
        if i<len(InputList)-1:
            OngoingString+=','
    OngoingString+=']'
    return OngoingString

def lowerNumberCount(InputList,InputNumber):
    # In a list of numbers, counts the occurrences of numbers lower than the current number
    i = 0
    for Number in InputList:
        if Number < InputNumber:
            i += 1
    return i

def makeBracketFigure(InputSequence,return_canvas=False,thioester='NHNH2',show_thioester_placeholder=True,show_special_aa_annotations=True,show_segment_labels=True,colors=(('#DF9282', '#926358'),('#EADA6C', '#A1964C'),('#AAE1DF', '#789D9C'),('#ADA9E4', '#77749C'),('#E7A4E0', '#9A748A'),('#E7A963', '#9F7646'),('#94D374', '#6C9A56'),('#AEC0E6', '#727E98'),('#C6A6DA', '#715D7F')),highlight_color='firebrick',show_rxn_label_placeholder=True,rxn_label_placeholder_text='NCL',aa_start_number=1,explicit_steps=True,good_thioesters=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ()-[]'),align_segments_at_top=True,highlight_poor_thioesters=True,px_per_aa=3,segment_spacing=20,desulf_pairs={"A":"C","V":"Pen"},highlight_non_cys_thiols=True):
    # This is the main function of BracketMaker. It takes input in the form of a correctly formatted text string, and returns either a drawSvg canvas object (return_canvas=True) or the contents of a SVG file (=False default)
    # Note- when I'm cleaning my code, a few of these variables are not being used anymore
    # Preview mode toggle; will change settings for multiple fxns
    global preview_mode
    preview_mode = return_canvas
    ## (All integer values below are in "points" for SVG or px for PNG; by default, these are equivalent
    # FLAG FIX - which of these need to be global variables, and which should be pushed into their sub-fxns?
    ## CANVAS DIMENSIONS
    global bottom_segment_width
   # bottom_segment_width = 300 # All segment sizes (and the entire canvas) will be calculated relative to this...moved to later
    canvas_margin = 10 # Blank space outside the border
    canvas_padding = 10 # Blank space inside the border
    ## CANVAS STYLE
    canvas_border_active = True
    border_stroke_width = 2
    border_stroke_color = '#000000'
    canvas_fill_color = '#ffffff'
    # SEQUENCE FIGURE
    sequence_font_size = 18
    show_sequence = True
    global sub_numbers
    sub_numbers = True # Numbers are subscripted in SVG file
    sequence_font = "Courier New" # NOTE: May need to change for PC vs. Mac
    padding_between_figures = 48
    global aa_start_number_global
    aa_start_number_global = aa_start_number # Have to assign to a second global variable for other fxns to use this
    # COLORS
    terminal_segment_stroke_color = "#444" # Can be same as default_segment_stroke_color, but it's useful to make terminal segments stand out
    global default_segment_fill_color
    default_segment_fill_color = "#D3D3D3" # Will be overwritten if show_segment_gradients = True
    global default_segment_stroke_color
    default_segment_stroke_color = "#808080" # Will be overwritten if show_stroke_colors = True
    global show_segment_gradients
    show_segment_gradients = True # If true, larger segments will have colorful gradients matching their internal pieces
    global show_stroke_colors
    show_stroke_colors = True # If true, this will use a specified list of stroke colors
    global default_font_color
    default_font_color = '#000000'
    # SEGMENT STYLE
    global explicit_steps_global
    explicit_steps_global = explicit_steps
    global highlight_color_global
    highlight_color_global = highlight_color
    global segment_stroke_width
    segment_stroke_width=2
    global segment_height
    segment_height = 30
    global segment_horizontal_spacing
    segment_horizontal_spacing = 30 # Adjustment factor; higher values will cause graph to be more "spread out" horizontally
    global segment_funnel_scaling
    segment_funnel_scaling=0.8 # Ratio for lower segments to give the diagram a nice funnel shape
    segment_label_font = "Arial" # Many standard fonts are supported.
    global segment_label_font_size
    segment_label_font_size = 16
    global annotation_font_size
    annotation_font_size = 12
    segment_corner_radius = 10 # Set to 0 for ordinary rectangles
    global segment_x_padding
    segment_x_padding = segment_spacing # Minimum white space between segments
    global px_per_aa_global
    px_per_aa_global=px_per_aa # Width of segments relative to AA length
    # RXN LINE STYLE
    global line_vertical_length
    line_vertical_length = 70
    line_stroke_width = 2
    line_stroke_color = "#808080"
    global rxn_label_font_size
    rxn_label_font_size = 12
    rxn_label_font = "Times New Roman" # Many standard fonts are supported.
    global rxn_label_padding
    rxn_label_padding = 10
    global rxn_label_placeholder_global
    rxn_label_placeholder_global = rxn_label_placeholder_text
    # ANNOTATIONS STYLE
    global annotation_line_color
    annotation_line_color="#AAAAAA"
    global annotation_stroke_width
    annotation_stroke_width=1.5

    # PREVIEW MODE SETTINGS; remove or simplify formatted text
    if return_canvas==True:
        show_sequence = False
        sub_numbers = False

    # SEQUENCE REFORMATTING
    # Check for valid characters; end fxn if invalid character is encountered
    format_tuple = checkFormat(InputSequence)
    if format_tuple[0] == False:
        print(f'Format error, type = {format_tuple[1]}')
        return
    # If correct format is detected, get list of segments and ligations from format checker
    global SegmentLigationList
    SegmentLigationList = format_tuple[2]
    # If there is a blank segment, this is actually a blank space between two explicit rxns; reformat for BM interpretation as a keyword with asterisk
    while '' in SegmentLigationList:
        i=SegmentLigationList.index('')
        FirstTag=stringToList(SegmentLigationList[i-1])
        SecondTag=stringToList(SegmentLigationList[i+1])
        ExplicitStepLabel = SecondTag[0]+'*'
        SecondTag[0]=ExplicitStepLabel
        CombinedTag=listToString(FirstTag+SecondTag)
        TempList = SegmentLigationList[0:i-1]+[CombinedTag]+SegmentLigationList[i+2:]
        SegmentLigationList = TempList
    # Update input sequence
    InputSequence=''.join(SegmentLigationList)

    # ANALYZE SEQUENCE
    # Split into lists of segments & ligation tags; the below calculations will reference these lists
    TerminalSegmentList = SegmentLigationList[0::2]
    LigationTagList = SegmentLigationList[1::2]

    # Make 'keyword match tracker' (for syntax checking; helps identify rxn steps that aren't matched to AAs)
    global KeywordMatchTracker
    KeywordMatchTracker = [] # Order corresponds to LigationTagList
    for i,tag in enumerate(LigationTagList):
        KeywordMatchTracker.append({})
        tag2 = stringToList(tag)
        if len(tag2)>1:
            tempdict = {}
            for term in tag2[1:]:
                term = term.rstrip('*')
                KeywordMatchTracker[i].update({term:False})

    #           LEVELS
    # List of all numbers in sequence
    global AllLevelNumbers
    AllLevelNumbers = re.findall(r"(\d+)", str(LigationTagList))
    # Turn this to a list of integers:
    TempList = []
    for Number in AllLevelNumbers:
        Number = int(Number)
        TempList.append(Number)
    AllLevelNumbers = TempList
    global BottomLevelNumber
    BottomLevelNumber = TopLevelNumber = 0
    if len(AllLevelNumbers)>0:
        BottomLevelNumber = min(AllLevelNumbers)
        TopLevelNumber = max(AllLevelNumbers)

    # C TERM LABELS
    global CTermLabelDict
    CTermLabelDict={}
    AACount = aa_start_number-1
    for Segment in TerminalSegmentList:
        # Determine if it has a C-terminal label; this is detected by a hyphen not inside of parentheses
        InsideParentheses=False
        AfterHyphen=False
        AfterHyphenList=[]
        for Char in Segment:
            if Char=='(':
                InsideParentheses=True
            if Char==')':
                InsideParentheses=False
            if AfterHyphen==True:
                AfterHyphenList.append(Char)
            if Char=='-' and InsideParentheses==False:
                AfterHyphen=True
            if InsideParentheses==False and AfterHyphen==False:
                AACount+=1
        # Add detected or default C-term label to dictionary; key is last AA number (NOT INDEX)
        # Don't add for final segment
        if len(AfterHyphenList)==0:
            CurrentThioester=thioester
            if show_thioester_placeholder==False:
                CurrentThioester=''
        else:
            CurrentThioester="".join(AfterHyphenList)
        CTermLabelDict.update({AACount:''})
        if CurrentThioester!="":
            CTermLabelDict.update({AACount:'-'+CurrentThioester})

    # For final segment, update dictionary to remove thioester if there isn't one specifically labeled; all variables should be set from last loop through
    if len(AfterHyphenList)==0:
        CTermLabelDict.update({AACount:''})
    # Finally, clean up our list of segments by removing thioesters (identified as hyphens not in parentheses)
    templist = []
    for Segment in TerminalSegmentList:
        CharList=[]
        AfterHyphen=InsideParentheses=False
        for Char in Segment:
            if Char=='(':
                InsideParentheses=True
            if Char==')' and InsideParentheses==True:
                InsideParentheses=False
            if Char=='-' and InsideParentheses==False:
                AfterHyphen=True
            if AfterHyphen==False:
                CharList.append(Char)
        StringMinusThioester = ''.join(CharList)
        templist.append(StringMinusThioester)
    TerminalSegmentList = templist

    # DEFINE SEQUENCE IN LIST FORM, PLUS SPECIAL AAS, PROTECTING GROUPS, ETC - each position contains info about each AA, including timing of special reactions
    # ALSO GET RANGES FOR OTHER FUNCTIONS
    global SequenceList
    SequenceList = []
    RecursiveRangeInput = '' # A string w/placeholders for characters and tags, to be called with main bracket function
    global SequenceRangeList
    SequenceRangeList = [] # Populate with ranges of each segment
    global ColorStopList
    ColorStopList = [] # Populate with a tuple specifying 'stop' point of each color/stroke pair in the list
    colors=colors*99
    c = 0 # Index of color list
    a = 0 # AA count (by index, for range list)
    # Based on desulfurization pairs dictionary, tell the function which AAs are affected by 'Desulfurization' keyword
    DesulfKeywordList=[]
    for after_key in desulf_pairs:
        before_key=desulf_pairs[after_key]
        SubString=f'{before_key}>{after_key}'
        DesulfKeywordList.append(SubString)
    # Begin looping through segments
    for iseg,Segment in enumerate(TerminalSegmentList):
        a0 = a # Left bound (index) of this segment; will match right bound of prev segment
        InsideParentheses=False
        for Char in Segment:
            if Char=='(':
                InsideParentheses=True
                InsideParenthesesList=[]
            # Normal AAs, simply append to sequence list
            if InsideParentheses==False:
                SequenceList.append(Char)
                RecursiveRangeInput+='X'
                a+=1
            # End of parentheses, identify special AA
            if Char==')' and InsideParentheses==True:
                InsideParentheses=False
                SpecialAAString = ''.join(InsideParenthesesList)
                # Special AA identified; determine what type it is
                listform = re.split('>',SpecialAAString) # Will be single length list unless explicit AA change is stated (e.g, C>A)
                CharIDDict = {}
                for n in range(BottomLevelNumber,TopLevelNumber+2): # Dictionary contains identity of this AA at all segment levels; initially assume it's equivalent to start AA
                    CharIDDict.update({n:listform[0]})
                pglistform = re.split('-',listform[0])
                # PROTECTING GROUPS - extract PG name and find level of removal step. Only applies to first term in list.
                if len(pglistform)>1:
                    PGName = '-'.join(pglistform[1:])
                    # Determine the level at which this PG is removed by searching adjacent tags
                    RemovalLevel = keywordSearch(PGName,LigationTagList,iseg)
                    # At removal level and lower, change to unprotected AA
                    for n in CharIDDict:
                        if n<=RemovalLevel:
                            CharIDDict[n]=pglistform[0]
                            listform[0]=pglistform[0] # Appropriately abbreviate keyword for the following step (useful for C-Acm>A)
                # AA CHANGES - find name of start and end AA and find level of change
                if len(listform)>1:
                    for iterm in range(1,len(listform)):
                        RxnName = listform[iterm-1]+'>'+listform[iterm]
                        RemovalLevel = keywordSearch(RxnName,LigationTagList,iseg)
                        # If this keyword is in our list of possible Desulfurization strings, also check for the Desulfurization keyword; use whichever is higher
                        if RxnName in DesulfKeywordList:
                            DesulfRemovalLevel = keywordSearch('Desulfurization',LigationTagList,iseg)
                            if DesulfRemovalLevel>RemovalLevel:
                                RemovalLevel=DesulfRemovalLevel
                        for n in CharIDDict:
                            if n<=RemovalLevel:
                                CharIDDict[n]=listform[iterm]
                # Define dictionary entry with its identity at all levels
                SequenceList.append(CharIDDict)
                RecursiveRangeInput+='X'
                a+=1
            # Inside parentheses - capture contents for keyword search later
            if InsideParentheses==True and Char!='(':
                InsideParenthesesList.append(Char)
        # End of segment (except last), add ligation tag to ongoing string
        if iseg<len(LigationTagList) and LigationTagList[iseg]!='[]':
            RecursiveRangeInput+=LigationTagList[iseg]
        # End of all segments, add range to range list
        SequenceRangeList.append((a0,a))
        # Add color of appropriate gradient to gradient list
        ColorStopList.append((a0,colors[c][0],colors[c][1]))
        c+=1

    # INITIALIZE OBJECT DICTIONARIES
    global SegmentDict
    SegmentDict = {}
    global SegmentIDCounter
    SegmentIDCounter = 1
    SequenceFigureDict = {}

    # Update bottom segment width based on number of segments:
    bottom_segment_width=100*len(TerminalSegmentList)
    if bottom_segment_width>500:
        bottom_segment_width=500

    # CALL MAIN BRACKET FUNCTION
    # Get bracket figure
    print("Getting bracket figure parameters...")
    drawBracket(recursiveRange(RecursiveRangeInput,BottomLevelNumber), BottomLevelNumber, good_thioesters=good_thioesters)
    print("Success")

    # BRACKET OBJECT ADJUSTMENTS & CANVAS DIMENSIONS
    canvas_adjustments = True # Turn False for development purposes only
    if canvas_adjustments == True:
        # ALIGN ALL TERMINAL SEGMENTS AT TOP
        if align_segments_at_top == True:
            # Check vertical position of all terminal segments, find highest one
            ymax=0
            for SegmentID in SegmentDict:
                if SegmentDict[SegmentID]['terminal']==True:
                    if SegmentDict[SegmentID]["y"]>ymax:
                        ymax=SegmentDict[SegmentID]["y"]
            # Check which need to be adjusted
            VerticalAlignList=[]
            for SegmentID in SegmentDict:
                if SegmentDict[SegmentID]['terminal']==True:
                    if SegmentDict[SegmentID]["y"]<ymax:
                        dy = ymax-SegmentDict[SegmentID]["y"]
                        VerticalAlignList.append((SegmentID,dy))
            # Make the necessary adjustments in dictionary
            for i,y in VerticalAlignList:
                oldy=SegmentDict[i]['y']
                oldlabely=SegmentDict[i]['labely']
                SegmentDict[i].update({'y':oldy+y,'labely':oldlabely+y})
        # INITIAL CANVAS VARIABLES
        canvas_x_left = canvas_x_right = canvas_y_top = 0
        # SHIFT OVERLAPPING SEGMENTS LEFT OR RIGHT
        # First, bin segments by height, also passing along center position which will be used for sorting
        TempList = [(n,SegmentDict[n]["center"],SegmentDict[n]["y"]) for n in SegmentDict]
        HeightBinDict = {}
        for (n,c,y) in TempList:
            if not y in HeightBinDict.keys():
                HeightBinDict.update({y:[]})
            HeightBinDict[y].append((n,c))
        # Create 'queue' with all of the layers we need to go down through
        LayerQueue=sorted(HeightBinDict.keys(),reverse=True)
        # LOOP THROUGH EACH LAYER FROM TOP DOWN
        TopLayer = LayerQueue[0]
        while len(LayerQueue)>0:
            # Remove the first entry from the queue, and get the corresponding segments at this layer
            CurrentHeight = LayerQueue[0]
            LayerQueue.remove(CurrentHeight)
            CurrentHeightSegments = HeightBinDict[CurrentHeight]
            # Find our center point. I've found that using the center line of canvas is not suitable; must anchor to the centermost segment.
            CentermostSegmentID = sorted(CurrentHeightSegments, key = lambda x: abs(x[1]))[0][0]
            # If this is the top layer, define our canvas dimensions based on this segment (assuming for now that it is the only object on the highest layer, only relevant for single segments really)
            if CurrentHeight==TopLayer:
                canvas_x_left=SegmentDict[CentermostSegmentID]["x"]
                canvas_x_right = SegmentDict[CentermostSegmentID]["x"] + SegmentDict[CentermostSegmentID]["width"] + SegmentDict[CentermostSegmentID]["ctermwidth"] + segment_x_padding
                canvas_y_top = SegmentDict[CentermostSegmentID]["y"] + SegmentDict[CentermostSegmentID]["height"] + 15
            # ID already sorts these left to right. Group into left list (to be moved further left) and right list (to be moved further right)
            # Don't move the center segment
            LeftIDList=[]
            RightIDList=[]
            for (n,c) in CurrentHeightSegments:
                if n<CentermostSegmentID:
                    LeftIDList.append(n)
                if n>CentermostSegmentID:
                    RightIDList.append(n)
            LeftIDList.sort(reverse=True)
            # Move segments in the current layer, and keep track of those we've moved
            MovedSegmentsRecord = []
            # Left segments - check if they clash, and move further left if so
            CompareID = CentermostSegmentID
            for myid in LeftIDList:
                # The right edge of this segment (including padding) should be less than or touching the left edge of the compared segment
                my_x_right = SegmentDict[myid]["x"] + SegmentDict[myid]["width"] + SegmentDict[myid]["ctermwidth"] + segment_x_padding
                your_x_left = SegmentDict[CompareID]["x"]
                if my_x_right>your_x_left:
                    x_shift = your_x_left-my_x_right # Will be negative
                    # Shift this segment and its labels over accordingly and update dictionary
                    oldx = SegmentDict[myid]["x"]
                    oldlabelx = SegmentDict[myid]["labelx"]
                    oldctermx = SegmentDict[myid]["ctermx"]
                    oldcenter = SegmentDict[myid]["center"] # checked by parent
                    newx = oldx + x_shift
                    newlabelx = oldlabelx + x_shift
                    newctermx = oldctermx + x_shift
                    newcenter = oldcenter + x_shift
                    SegmentDict[myid].update({"x":newx})
                    SegmentDict[myid].update({"labelx":newlabelx})
                    SegmentDict[myid].update({"ctermx":newctermx})
                    SegmentDict[myid].update({"center":newcenter})
                    # Record that we have moved this segment, and need to check the attached objects
                    my_parent_id = SegmentDict[myid]["parent"]
                    MovedSegmentsRecord.append((myid,my_parent_id,x_shift))
                # If this is further left than the current canvas dimensions, adjust boundaries accordingly
                my_x_left = SegmentDict[myid]["x"]
                if my_x_left<canvas_x_left:
                    canvas_x_left = my_x_left
                # Compare next segment to this one
                CompareID = myid
            # Right segments - determine if they clash, and shift further right
            CompareID = CentermostSegmentID
            for myid in RightIDList:
                # The left edge of this segment (including padding) should be touching or greater than the right edge of the compared segment
                my_x_left = SegmentDict[myid]["x"]
                your_x_right = SegmentDict[CompareID]["x"] + SegmentDict[CompareID]["width"] + SegmentDict[CompareID]["ctermwidth"] + segment_x_padding
                if my_x_left<your_x_right:
                    x_shift = your_x_right-my_x_left # Will be positive
                    # Shift this segment and its labels over accordingly and update dictionary
                    oldx = SegmentDict[myid]["x"]
                    oldlabelx = SegmentDict[myid]["labelx"]
                    oldctermx = SegmentDict[myid]["ctermx"]
                    oldcenter = SegmentDict[myid]["center"]
                    newx = oldx + x_shift
                    newlabelx = oldlabelx + x_shift
                    newctermx = oldctermx + x_shift
                    newcenter = oldcenter + x_shift
                    SegmentDict[myid].update({"x":newx})
                    SegmentDict[myid].update({"labelx":newlabelx})
                    SegmentDict[myid].update({"ctermx":newctermx})
                    SegmentDict[myid].update({"center":newcenter})
                # If this is further right than the current canvas dimensions, adjust boundaries accordingly
                my_x_right = SegmentDict[myid]["x"] + SegmentDict[myid]["width"] + SegmentDict[myid]["ctermwidth"] + segment_x_padding
                if my_x_right>canvas_x_right:
                    canvas_x_right = my_x_right
                # Compare next segment to this one
                CompareID = myid
            # Segments with children - make sure parent segment & rxn line are centered underneath daughter segments
            for (myid,c) in CurrentHeightSegments:
                if 'leftid' in SegmentDict[myid]:
                    # Get L and R bounds of Left and Right daughter segments
                    LeftSegmentDict=SegmentDict[SegmentDict[myid]["leftid"]]
                    LeftBound=LeftSegmentDict["x"]
                    RightSegmentDict=SegmentDict[SegmentDict[myid]["rightid"]]
                    RightBound=RightSegmentDict["x"]+RightSegmentDict["width"]+RightSegmentDict["ctermwidth"]
                    # Center bracket arms underneath daughter segments
                    LeftCenter = LeftSegmentDict["center"]
                    MyLeftLine = SegmentDict[myid]["linexl"]
                    if MyLeftLine!=LeftCenter:
                        SegmentDict[myid].update({"linexl":LeftCenter})
                    RightCenter = RightSegmentDict["center"]
                    MyRightLine = SegmentDict[myid]["linexr"]
                    if MyRightLine!=RightCenter:
                        SegmentDict[myid].update({"linexr":RightCenter})
                    # Center parent segment and rxn line underneath daughter segments
                    MidBound=(LeftBound+RightBound)/2
                    MyCenter=SegmentDict[myid]["center"]
                    x_shift=MidBound-MyCenter
                    if x_shift!=0:
                        oldx = SegmentDict[myid]["x"]
                        oldlabelx = SegmentDict[myid]["labelx"]
                        oldctermx = SegmentDict[myid]["ctermx"]
                        oldcenter = SegmentDict[myid]["center"]
                        newx = oldx + x_shift
                        newlabelx = oldlabelx + x_shift
                        newctermx = oldctermx + x_shift
                        newlinexc = newx+(SegmentDict[myid]["width"]/2)
                        newrxnlabelx = newlinexc+rxn_label_padding
                        newcenter = oldcenter + x_shift
                        SegmentDict[myid].update({"x":newx})
                        SegmentDict[myid].update({"labelx":newlabelx})
                        SegmentDict[myid].update({"ctermx":newctermx})
                        SegmentDict[myid].update({"linexc":newlinexc})
                        SegmentDict[myid].update({"rxnlabelx":newrxnlabelx})
                        SegmentDict[myid].update({"center":newcenter})
                elif 'centerid' in SegmentDict[myid]:
                    # Get L and R bounds of daughter segment
                    CenterSegmentDict=SegmentDict[SegmentDict[myid]["centerid"]]
                    YourCenter=CenterSegmentDict["center"]
                    MyCenter=SegmentDict[myid]["center"]
                    x_shift=YourCenter-MyCenter
                    # Shift this segment and its label over accordingly
                    if x_shift!=0:
                        oldx = SegmentDict[myid]["x"]
                        oldlabelx = SegmentDict[myid]["labelx"]
                        oldctermx = SegmentDict[myid]["ctermx"]
                        oldcenter = SegmentDict[myid]["center"]
                        newx = oldx + x_shift
                        newlabelx = oldlabelx + x_shift
                        newctermx = oldctermx + x_shift
                        newsinglelinex = YourCenter
                        newrxnlabelx = newsinglelinex + rxn_label_padding
                        newcenter = oldcenter + x_shift
                        SegmentDict[myid].update({"x":newx})
                        SegmentDict[myid].update({"labelx":newlabelx})
                        SegmentDict[myid].update({"ctermx":newctermx})
                        SegmentDict[myid].update({"singlelinex":newsinglelinex})
                        SegmentDict[myid].update({"rxnlabelx":newrxnlabelx})
                        SegmentDict[myid].update({"center":newcenter})

    # CALCULATE CANVAS DIMENSIONS
    canvas_height = canvas_y_top
    canvas_width = canvas_x_right-canvas_x_left
    canvas_x_center = (canvas_x_right+canvas_x_left)/2

    # SEQUENCE FIGURE
    # Get sequence figure (and adjust top of canvas to match)
    if show_sequence == True:
        print("Getting sequence figure parameters...")
        # Initialize variables
        AACount = 0
        SequenceTspanList = []
        UnusedAAList = list('123456789')
        EncounteredSpecialAA = {} # Dictionary of special AAs in legend
        UnusedAAIndex = 0
        SequenceDeltaY = sequence_font_size*2
        IsFirstTime = True
        #InsideParentheses = False
        # IsUnderlined = False
        LastColor='' # Trigger to start new color tspan
        CloseHighlight=False # Trigger to close highlight tspan
        LegendString=''
        # SEQUENCE FIGURE INITIAL CALCULATIONS
        # Width of the text block and position of the sequence will all be dependent on coordinates from drawBracket
        CharacterWidth = sequence_font_size*2/3 # (I think this is the ratio for Courier)
        CharacterLimit = int(canvas_width/CharacterWidth)
        # Round down to a multiple of 10
        SequenceRowLength = CharacterLimit - (CharacterLimit % 10)
        if SequenceRowLength==0:
            SequenceRowLength=10
        # Calculate number of rows total
        if len(SequenceList) % SequenceRowLength == 0:
            NumberOfRows = int(len(SequenceList)/SequenceRowLength)
        else:
            NumberOfRows = int(len(SequenceList)/SequenceRowLength+1)
        # Calculate height and adjust canvas top
        sequence_figure_padding=15
        SequenceHeight = sequence_font_size*2*NumberOfRows
        canvas_height+=(SequenceHeight+sequence_figure_padding)
        # All coordinates of objects within this figure are relative to the top left point, set here to 0,0.
        sequence_x_position=numbering_x_position=canvas_x_left
        numbering_y_position=canvas_height
        sequence_y_position=numbering_y_position-sequence_font_size
        # At beginning, begin new tspan tag with no vertical shift:
        SequenceTspanList.append(f'<tspan x="{sequence_x_position}" dy="0">')
        # DEFINE COLORS AT EACH POSITION with new dictionary
        CharacterColorDict={}
        for i,(aa,fill,stroke) in enumerate(ColorStopList):
            RightBound = len(SequenceList)
            if i+1<len(ColorStopList):
                RightBound=ColorStopList[i+1][0]
            for n in range(aa,RightBound):
                CharacterColorDict.update({n:fill})
        # LOOP THROUGH CHARACTERS TO DRAW SEQUENCE FIGURE
        for i,aa in enumerate(SequenceList):
            # Determine if the start of this segment is also the start of a row
            IsRowStart = False
            if AACount % SequenceRowLength == 0:
                IsRowStart = True
            # Add tspan tag for new color
            if CharacterColorDict[i]!=LastColor:
                # Except for first time, add ending tspan
                if IsFirstTime == False:
                    SequenceTspanList.append(r'</tspan>')
                IsFirstTime=False
                SequenceTspanList.append(f'<tspan style="fill: {CharacterColorDict[i]};">')
                # Remember this color to check next letter
                LastColor = CharacterColorDict[i]
            # Add characters with appropriate spacing and highlighting
            # If this is a special AA, decide what to do with it
            if isinstance(aa,dict):
                # Show whatever is in the final sequence
                aa=aa[BottomLevelNumber]
                # If more than 1 character, assign an abbreviation and record in legend
                if len(aa)>1:
                    aa_to_check=aa
                    if aa_to_check in EncounteredSpecialAA: # If encountered before, use the existing dict key; legend will already have this
                        aa=EncounteredSpecialAA[aa_to_check]
                    else:
                        if len(LegendString)>0:
                            LegendString+=', '
                        LegendString+=f'{UnusedAAList[UnusedAAIndex]} = {aa}'
                        aa=UnusedAAList[UnusedAAIndex]
                        UnusedAAIndex+=1
                        EncounteredSpecialAA[aa_to_check]=aa # update dictionary for later
                # Add tspan tag for highlight
                SequenceTspanList.append(f'<tspan text-decoration="underline" fill="{highlight_color}">')
                CloseHighlight=True
            # Add AA
            SequenceTspanList.append(aa)
            if CloseHighlight==True:
                SequenceTspanList.append('</tspan>')
                CloseHighlight=False
            # Add a space after every tenth character
            if (i+1) % 10 == 0:
                SequenceTspanList.append(" ")
            # Close tspan tags and open new ones after every row
            if (i+1) % SequenceRowLength == 0:
                # End previous tspan tags (color and position)
                SequenceTspanList.append(r'</tspan></tspan>')
                # Add next tspan tags (color and position)
                SequenceTspanList.append(f'<tspan x="{sequence_x_position}" dy="{SequenceDeltaY}"><tspan style="fill: {CharacterColorDict[i]};">')
        # Add legend and closing tags
        SequenceTspanList.append('</tspan>') # Closes color tag; text following this will be black
        if len(LegendString)>0:
            # Calculate whether this will be placed in the current or next line based on remaining white space
            RemainingSpaces = SequenceRowLength - ((i+1) % SequenceRowLength)
            # Include white spaces that would have been placed every 10th character
            RemainingSpaces = RemainingSpaces + int(RemainingSpaces/10)
            if len(LegendString)<=RemainingSpaces:
                LegendSpacing = RemainingSpaces - len(LegendString)
                SequenceTspanList.append(r'<tspan font-weight="normal">')
                SequenceTspanList.append("&#160;"*LegendSpacing+LegendString)
                SequenceTspanList.append(r'</tspan>')
            # Otherwise add to the next line, and adjust the position of the entire figure and canvas
            else:
                SequenceTspanList.append(f'<tspan x="{sequence_x_position}" dy="{SequenceDeltaY}" font-weight="normal">')
                SequenceTspanList.append(LegendString)
                SequenceTspanList.append(r'</tspan>')
                canvas_height+=SequenceDeltaY
                sequence_y_position+=SequenceDeltaY
        # Add final ending tspan for line height
        SequenceTspanList.append(r'</tspan>')
        # FINISH AND ADD TO FIGURE
        # Condense to text string
        FormattedSequence = ''.join(SequenceTspanList)
        # NUMBERING - this adds numbers above sequence, labeling above every 10th character
        # Initialize variables
        NumberingTspanList = []
        GlobalNumCount = 0
        # Add starting tspan with no vertical shift
        NumberingTspanList.append(f'<tspan x="{sequence_x_position}" dy="0">')
        # Determine starting and ending number
        NumberingStart = aa_start_number + 9
        NumberingEnd = aa_start_number - 1 + (len(SequenceList) - (len(SequenceList) % 10))
        # Loop through all numbers in this range
        for Number in range(NumberingStart,NumberingEnd+1,10):
            DigitCount = len(str(Number))
            # Add an appropriate number of spaces, then the number
            Spaces = '&#160;'*(10-DigitCount)
            NumberingTspanList.append(Spaces)
            NumberingTspanList.append(str(Number)+'&#160;')
            # Update our counter
            GlobalNumCount += 10
            # At end of row, close tspan tag and start next, and reset counter
            if GlobalNumCount >= SequenceRowLength:
                NumberingTspanList.append(f'</tspan>')
                GlobalNumCount = 0
                NumberingTspanList.append(f'<tspan x="{sequence_x_position}" dy="{SequenceDeltaY}">')
        # Add final ending tspan
        NumberingTspanList.append(r'</tspan>')
        # Condense to text string
        FormattedNumbering = ''.join(NumberingTspanList)
        # UPDATE OBJECT DEFINITION
        # Store all variables in global sequence figure dictionary
        SequenceFigureDict.update({
            "seq":FormattedSequence,
            "seqx":sequence_x_position,
            "seqy":sequence_y_position,
            "num":FormattedNumbering,
            "numx":numbering_x_position,
            "numy":numbering_y_position
        })

    # DEFINE CANVAS
    total_document_width = canvas_width+(canvas_padding*2)
    total_document_height = canvas_height+(canvas_padding*3)
    xorigin = -(total_document_width/2-canvas_x_center)
    yorigin = -canvas_padding
    canvas = draw.Drawing(total_document_width, total_document_height, origin=(xorigin,yorigin), displayInline=False)
    # For development - draw a canvas border
    canvas_border = False
    if canvas_border == True:
        canvas.append(
            draw.Rectangle(canvas_x_left,0,canvas_width,canvas_height,fill=canvas_fill_color,stroke_width=2,stroke='#000000')
                )

    # SEQUENCE FIGURE AND NUMBERING
    if show_sequence == True:
        SequenceFigure = draw.Text(
            SequenceFigureDict["seq"],
            sequence_font_size,
            SequenceFigureDict["seqx"],
            SequenceFigureDict["seqy"],
            font_family=sequence_font,
            font_weight="bold")
        NumberingFigure = draw.Text(
            SequenceFigureDict["num"],
            sequence_font_size,
            SequenceFigureDict["numx"],
            SequenceFigureDict["numy"],
            font_family=sequence_font,
            font_weight="bold",
            fill=default_font_color)
        canvas.append(SequenceFigure)
        canvas.append(NumberingFigure)

    # SEGMENTS AND RXN LINES; loop through Segment dictionary
    for i in SegmentDict:
        # Internal segments, define gradients for fill and/or stroke:
        FillColor = SegmentDict[i]["fill"]
        if isinstance(FillColor, list):
            # Define gradient vector (x1, y1, x2, y2), a line stretching from one end of segment to the other. Y doesn't matter as long as it's the same.
            FillGradient = draw.LinearGradient(SegmentDict[i]["x"], SegmentDict[i]["y"], SegmentDict[i]["x"]+SegmentDict[i]["width"], SegmentDict[i]["y"])
            # Add gradient stops as defined in list
            prevcolor = ''
            for stop,color in FillColor:
                if prevcolor!='':
                    FillGradient.addStop(stop,prevcolor)
                FillGradient.addStop(stop,color)
                prevcolor=color
            # Update fill color to be this drawSVG gradient object
            SegmentDict[i]["fill"] = FillGradient
        StrokeColor = SegmentDict[i]["stroke"]
        if isinstance(StrokeColor, list):
            # Define gradient vector (x1, y1, x2, y2), a line stretching from one end of segment to the other. Y doesn't matter as long as it's the same.
            StrokeGradient = draw.LinearGradient(SegmentDict[i]["x"], SegmentDict[i]["y"], SegmentDict[i]["x"]+SegmentDict[i]["width"], SegmentDict[i]["y"])
            # Add gradient stops as defined in list
            prevcolor = ''
            for stop,color in StrokeColor:
                if prevcolor!='':
                    StrokeGradient.addStop(stop,prevcolor)
                StrokeGradient.addStop(stop,color)
                prevcolor=color
            # Update stroke color to be this drawSVG gradient object
            SegmentDict[i]["stroke"] = StrokeGradient
        # Define shape and add to canvas
        segment_shape = draw.Rectangle(
            SegmentDict[i]["x"],
            SegmentDict[i]["y"],
            SegmentDict[i]["width"],
            SegmentDict[i]["height"],
            fill=SegmentDict[i]["fill"],
            stroke_width=SegmentDict[i]["strokewidth"],
            stroke = SegmentDict[i]["stroke"],
            rx=segment_corner_radius,
            fill_opacity=SegmentDict[i]["opacity"])
        canvas.append(segment_shape)

        # SEGMENT LABELS
        # At this point, the function has information about special AAs that need to be included on each segment
        # Central label will contain first and last AA spaced out with an appropriate amount of dots
        # First, calculate the number of dots - enough to push first and last AA to the edges of segment
        SegmentWidth=SegmentDict[i]["width"] # Width of our bounding box
        StringToCheck=f'{SegmentDict[i]["firstAA"]}{SegmentDict[i]["firstAAnumber"]}{SegmentDict[i]["lastAA"]}{SegmentDict[i]["lastAAnumber"]}' # e.g., C124R158
        WidthOfLabelText=getCharWidth(StringToCheck,FontSize=segment_label_font_size,sub_numbers=not preview_mode) # Width of text we need to space out
        SpaceToFill=SegmentWidth-WidthOfLabelText
        # Terminal segments have an inner stroke, so the inside space is just slightly less
        if SegmentDict[i]['terminal']==True:
            SpaceToFill-=3 # Should be same as terminal segment stroke width
        dotwidth=0.28*segment_label_font_size # Found experimentally, measuring in canvas
        DotsToAdd=int(SpaceToFill/dotwidth) # Will correctly round up/down, but still save as integer. Sometimes that extra 1 dot makes a difference.
        if DotsToAdd<3:
            DotsToAdd=3 # If there are less than 3 dots worth of space, label will look too large for segment; fix by using a smaller font size, or larger segment width
        spacer='.'*DotsToAdd # Define spacer
        # Now, begin building the string for the segment label
        # FIRST AA
        segment_label_string=""
        if SegmentDict[i]["highlightfirstAA"]==True and preview_mode==False:
            segment_label_string+=f'<tspan fill = "{highlight_color}">' # Optional opening tspan for highlight
        segment_label_string+=SegmentDict[i]["firstAA"]+str(SegmentDict[i]["firstAAnumber"]) # First AA and number
        if SegmentDict[i]["highlightfirstAA"]==True and preview_mode==False:
            segment_label_string+=f'</tspan>' # Optional closing tspan for highlight
        # SPACER
        segment_label_string+=spacer # Spacer in between AAs, calculated above
        # LAST AA
        if SegmentDict[i]["highlightlastAA"]==True and preview_mode==False:
            segment_label_string+=f'<tspan fill = "{highlight_color}">' # Optional opening tspan for highlight
        segment_label_string+=SegmentDict[i]["lastAA"]+str(SegmentDict[i]["lastAAnumber"]) # Last AA and number
        if SegmentDict[i]["highlightlastAA"]==True and preview_mode==False:
            segment_label_string+=f'</tspan>' # Optional closing tspan for highlight

        # In SVG, subscript numbers
        if preview_mode==False:
            segment_label_string=subNumbers(segment_label_string)
        segment_label = draw.Text(
            segment_label_string,
            segment_label_font_size,
            SegmentDict[i]["labelx"],
            SegmentDict[i]["labely"],
            center=True,
            valign="bottom",
            lineHeight=1,
            font_family=segment_label_font,
            font_weight="bold",
            fill=default_font_color)
        if show_segment_labels == True:
            canvas.append(segment_label)


        # ANNOTATIONS ABOVE SEGMENTS
        annotation_line_length=14
        segment_label_simplified=f'{SegmentDict[i]["firstAA"]}{SegmentDict[i]["firstAAnumber"]}{spacer}{SegmentDict[i]["lastAA"]}{SegmentDict[i]["lastAAnumber"]}' # Does not include tspan tags
        # Get the y-coordinate of top of segment, and width and position of label
        AnnotationStart=SegmentDict[i]["y"]+SegmentDict[i]["height"]-6
        LabelCenter=SegmentDict[i]["labelx"]
        LabelWidth=getCharWidth(segment_label_simplified,FontSize=segment_label_font_size,sub_numbers=not preview_mode)
        # Get the position of first AA label text (centered over AA + number)
        LabelLeftBorder=LabelCenter-(LabelWidth/2)
        FirstAALabelWidth=getCharWidth(f'{SegmentDict[i]["firstAA"]}',FontSize=segment_label_font_size,sub_numbers=not preview_mode)
        FirstAAPosition=LabelLeftBorder+FirstAALabelWidth/2
        # Get the position of the last AA label text (centered over AA + number)
        LabelRightBorder=LabelCenter+(LabelWidth/2)
        LastAALabelWidth=getCharWidth(f'{SegmentDict[i]["lastAA"]}',FontSize=segment_label_font_size,sub_numbers=not preview_mode)
        LastAAPosition=LabelRightBorder-LastAALabelWidth/2

        # FIRST AND LAST AA - Determine if this is a protected AA, and if so, add the appropriate annotation
        AnnotationDict=SegmentDict[i]["annotations"]
        FirstAANumber=SegmentDict[i]["firstAAnumber"]
        LastAANumber=SegmentDict[i]["lastAAnumber"]
        if FirstAANumber in AnnotationDict:
            # Draw line above first AA
            first_aa_line = draw.Lines(
                FirstAAPosition,
                AnnotationStart,
                FirstAAPosition,
                AnnotationStart+annotation_line_length,
                stroke_width=annotation_stroke_width, # Defined in user input section RXN LINE STYLE
                stroke = SegmentDict[i]["ntermstroke"], # Defined in user input section RXN LINE STYLE
                fill_opacity = 0, # No fill
                close=False)
            canvas.append(first_aa_line)
            # Add text box with protecting group
            FirstAnnotationText=AnnotationDict[FirstAANumber]
            first_aa_annotation = draw.Text(
                FirstAnnotationText,
                annotation_font_size,
                FirstAAPosition,
                AnnotationStart+annotation_line_length+2,
                center=True,
                valign="bottom",
                lineHeight=1,
                font_family=segment_label_font,
                font_weight="bold",
                fill=SegmentDict[i]["ntermfill"])
            if show_special_aa_annotations==True:
                canvas.append(first_aa_annotation)
        if LastAANumber in AnnotationDict:
            # Draw line above last AA
            last_aa_line = draw.Lines(
                LastAAPosition,
                AnnotationStart,
                LastAAPosition,
                AnnotationStart+annotation_line_length,
                stroke_width=annotation_stroke_width, # Defined in user input section RXN LINE STYLE
                stroke = SegmentDict[i]["ctermstroke"], # Defined in user input section RXN LINE STYLE
                fill_opacity = 0, # No fill
                close=False)
            canvas.append(last_aa_line)
            # Add text box with protecting group
            LastAnnotationText=AnnotationDict[LastAANumber]
            last_aa_annotation = draw.Text(
                LastAnnotationText,
                annotation_font_size,
                LastAAPosition,
                AnnotationStart+annotation_line_length+2,
                center=True,
                valign="bottom",
                lineHeight=1,
                font_family=segment_label_font,
                font_weight="bold",
                fill=SegmentDict[i]["ctermfill"])
            if show_special_aa_annotations==True:
                canvas.append(last_aa_annotation)
        # MIDDLE AA ANNOTATIONS - add annotation, but adjust to avoid running into other canvas objects
        AnnotationLeftBound=SegmentDict[i]["x"]
        AnnotationRightBound=SegmentDict[i]["x"]+SegmentDict[i]["width"]
        # Save this width for calculating relative line position later
        AnnotationBoundWidth=AnnotationRightBound-AnnotationLeftBound
        # Get a list of middle AAs that need to be annotated; these will be processed left-to-right
        MiddleTagList=[] # Contains (AA, number, lineposition) for each annotated AA
        for AANumber in AnnotationDict.keys():
            if AANumber!=FirstAANumber and AANumber!=LastAANumber:
                MiddleAnnotationText=AnnotationDict[AANumber]+str(AANumber)
                # Get relative position of this middle AA within acceptable boundaries
                MiddleLinePosition = AnnotationLeftBound+(AnnotationBoundWidth*(AANumber-FirstAANumber)/(LastAANumber-FirstAANumber))
                # Keep track of this line position and text box contents; tags will be combined later if necessary
                MiddleTagList.append((AANumber, MiddleAnnotationText, MiddleLinePosition))
        # Loop through list of tags, and determine if they should be placed as-is, or if they should be combined with the previous tag to avoid clashing
        TagsWithOwnBoxes=MiddleTagList[:] # List of tags that will get their own unique text box label
        RxnLineXShift=0 # For later, record if any of these text boxes clash with the central rxn line
        for (AANumber, MiddleAnnotationText, MiddleLinePosition) in MiddleTagList:
            TagTuple=(AANumber, MiddleAnnotationText, MiddleLinePosition)
            if TagTuple in TagsWithOwnBoxes:
                # If in the list, we are starting a new annotation tag, and checking to see which are combined with this
                TagsWithOwnBoxes.remove(TagTuple) # Remove tag so we don't check it against itself
                CurrentAnnotationList=[(TagTuple)] # Begin list of all tags that will be combined into this annotation, starting with the current tag
                # All middle tags will be shifted up relative to the N- and C-terminal tags, to avoid having to calculate for that clash
                ShiftUp=False
                if FirstAANumber in AnnotationDict or LastAANumber in AnnotationDict:
                    ShiftUp=True
                # Loop through tags which potentially have their own boxes, and determine if any need to be absorbed into this tag instead
                MyRightBound=MiddleLinePosition+(getCharWidth(MiddleAnnotationText,sub_numbers=not preview_mode))/2
                TagsToRemove=[] # Cannot remove items from TagsWithOwnBoxes during the below loop, so I have to do it after
                for z,(YourAANumber, YourText, YourLinePosition) in enumerate(TagsWithOwnBoxes):
                    ClashingTag=(YourAANumber, YourText, YourLinePosition)
                    Clash=False
                    # Determine if there is a clash between my right bound and your left bound
                    YourLeftBound=YourLinePosition-(getCharWidth(YourText,sub_numbers=not preview_mode))/2
                    if YourLeftBound<=MyRightBound:
                        Clash=True
                    # Determine if there will be a clash if all the next tags were to be absorbed into one
                    if Clash==False:
                        CombinedNextTags=condenseAnnotation(TagsWithOwnBoxes[z:])
                        PreviewCenterPosition=CombinedNextTags[2]
                        PreviewCombinedText=CombinedNextTags[1]
                        PreviewLeftBound=PreviewCenterPosition-getCharWidth(PreviewCombinedText,sub_numbers=not preview_mode)/2
                        if PreviewLeftBound<=MyRightBound:
                            Clash=True
                    # If so, the offending tag gets absorbed instead of getting its own box
                    if Clash==True:
                        # TagsWithOwnBoxes.remove(ClashingTag) # Can't do this while actively looping through the list
                        TagsToRemove.append(ClashingTag)
                        CurrentAnnotationList.append(ClashingTag)
                        # For the next loop checking the next possible clash, update our R boundary by calculating a new center position and text box width, using a custom small fxn
                        CondensedTag=condenseAnnotation(CurrentAnnotationList)
                        CondensedTagText=CondensedTag[1]
                        AvgLinePosition=CondensedTag[2]
                        MyRightBound=AvgLinePosition+(getCharWidth(CondensedTagText,sub_numbers=not preview_mode))/2
                # For all the tags we absorbed into the current, remove them from the list so they don't get their own box
                for x in TagsToRemove:
                    TagsWithOwnBoxes.remove(x)
                # After that loop, we now have a list of all the lines that need to be drawn and all the tag labels that need to be combined above it
                # First, draw all the lines; calculate line length based on ShiftUp toggle
                ThisLineLength=annotation_line_length
                if ShiftUp==True:
                    ThisLineLength=annotation_line_length*1.5
                # Get a little dictionary with the stroke color of each line, matching the appropriate stroke color of its segment
                AnnotationStrokeDict={}
                for Annotation in CurrentAnnotationList:
                    AAIndex=Annotation[0]-1
                    for (colorindex,fill,stroke) in ColorStopList:
                        if colorindex<=AAIndex:
                            AnnotationStrokeDict.update({Annotation[2]:stroke})
                # Define each line and add to canvas
                TheseLines=[x[2] for x in CurrentAnnotationList]
                for MiddleLinePosition in TheseLines:
                    # Change line color to be equal to corresponding stroke color
                    ThisLineColor=AnnotationStrokeDict[MiddleLinePosition]
                    middle_aa_line = draw.Lines(
                        MiddleLinePosition,
                        AnnotationStart,
                        MiddleLinePosition,
                        AnnotationStart+ThisLineLength,
                        stroke_width=annotation_stroke_width, # Defined in user input section RXN LINE STYLE
                        stroke=ThisLineColor, # Defined in user input section RXN LINE STYLE
                        fill_opacity = 0, # No fill
                        close=False)
                    canvas.append(middle_aa_line)

                # Get average text box position and condensed label
                CondensedTag=condenseAnnotation(CurrentAnnotationList)
                CondensedTagText=CondensedTag[1]
                CondensedTagWidth=getCharWidth(CondensedTagText, FontSize=annotation_font_size, sub_numbers=not preview_mode)
                if preview_mode==False:
                    CondensedTagText=subNumbers(CondensedTagText,FontSize=annotation_font_size)
                AvgLinePosition=CondensedTag[2]
                # Record for later - adjust rxn line horizontally so it does not clash with tags across the center
                # First, determine if this annotation crosses the segment center
                SegmentXCenter = SegmentDict[i]["x"]+(SegmentDict[i]["width"]/2)
                MyLeftBound = AvgLinePosition - (CondensedTagWidth/2)
                MyRightBound = AvgLinePosition + (CondensedTagWidth/2)
                leftxthreshold = SegmentDict[i]["x"] + 10 # Make sure rxn line doesn't spill off edge
                rightxthreshold = SegmentDict[i]["x"] + SegmentDict[i]["width"] - 10
                if MyLeftBound<=SegmentXCenter<=MyRightBound:
                    # Determine if this should shift left or right, depending which edge is closer
                    leftdistance = SegmentXCenter - MyLeftBound
                    rightdistance = MyRightBound - SegmentXCenter
                    if leftdistance<rightdistance and SegmentXCenter-leftdistance>=leftxthreshold:
                        RxnLineXShift = -leftdistance
                    elif rightdistance<=leftdistance and SegmentXCenter+rightdistance<=rightxthreshold:
                        RxnLineXShift = rightdistance
                # Calculate text box height based on ShiftUp toggle
                ThisBoxHeight=annotation_line_length+5
                if ShiftUp==True:
                    ThisBoxHeight=annotation_line_length*1.5+5
                if preview_mode==False:
                    # Add extra space for subscript numbers
                    ThisBoxHeight+=4
                # Define and add this text box
                middle_aa_annotation = draw.Text(
                    CondensedTagText,
                    annotation_font_size,
                    AvgLinePosition,
                    AnnotationStart+ThisBoxHeight,
                    center=True,
                    valign="bottom",
                    lineHeight=1,
                    font_family=segment_label_font,
                    font_weight="bold",
                    fill=default_font_color)
                if show_special_aa_annotations==True:
                    canvas.append(middle_aa_annotation)

        # Define C-term thioester label and add to canvas if it exists
        if SegmentDict[i]["cterm"] != "":
            cterm_label = draw.Text(
                SegmentDict[i]["cterm"],
                annotation_font_size,
                SegmentDict[i]["ctermx"],
                SegmentDict[i]["labely"],
                font_family=segment_label_font,
                font_weight="bold",
                fill=SegmentDict[i]["ctermfill"])
            canvas.append(cterm_label)
        # Define rxn lines and add to canvas if this segment has one defined
        if "linexl" in SegmentDict[i]:
            # First check if this is a bad thioester; if so, temporarily change the color we use to highlight
            usual_line_color=line_stroke_color
            usual_line_width=line_stroke_width
            if SegmentDict[i]["goodthioester"] == False and highlight_poor_thioesters==True:
                line_stroke_color=highlight_color
                line_stroke_width=usual_line_width*2
            if SegmentDict[i]["cysthiol"]==False and highlight_non_cys_thiols==True:
                line_stroke_color="#0008ff"
                line_stroke_width=usual_line_width*2
            # Adjust rxn line X shift to accommodate width of stroke
            if RxnLineXShift>0:
                RxnLineXShift+=line_stroke_width
            elif RxnLineXShift<0:
                RxnLineXShift-=line_stroke_width
            # Draw rxn lines as defined
            rxn_line_1 = draw.Lines(
                SegmentDict[i]["linexc"]+RxnLineXShift, SegmentDict[i]["lineyb"],
                SegmentDict[i]["linexc"]+RxnLineXShift, SegmentDict[i]["lineyc"],
                stroke_width=line_stroke_width, # Defined in user input section RXN LINE STYLE
                stroke = line_stroke_color, # Defined in user input section RXN LINE STYLE
                fill_opacity = 0, # No fill
                close=False)
            rxn_line_2 = draw.Lines(
                SegmentDict[i]["linexl"], SegmentDict[i]["lineyt"],
                SegmentDict[i]["linexl"], SegmentDict[i]["lineyc"],
                SegmentDict[i]["linexr"], SegmentDict[i]["lineyc"],
                SegmentDict[i]["linexr"], SegmentDict[i]["lineyt"],
                stroke_width=line_stroke_width, # Defined in user input section RXN LINE STYLE
                stroke = line_stroke_color, # Defined in user input section RXN LINE STYLE
                fill_opacity = 0, # No fill
                close=False)
            canvas.append(rxn_line_1)
            canvas.append(rxn_line_2)
            # Add rxn line extension for vertically-aligned terminal segments at the top
            if align_segments_at_top==True:
                for x,y in VerticalAlignList: # SegmentID, deltaY
                    if SegmentDict[i]["leftid"] == x:
                        rxn_line_3 = draw.Lines(
                            SegmentDict[i]["linexl"], SegmentDict[i]["lineyt"],
                            SegmentDict[i]["linexl"], SegmentDict[i]["lineyt"]+y,
                            stroke_width=line_stroke_width, # Defined in user input section RXN LINE STYLE
                            stroke = line_stroke_color, # Defined in user input section RXN LINE STYLE
                            stroke_dasharray = "10,10",
                            fill_opacity = 0, # No fill
                            close=False)
                        canvas.append(rxn_line_3)
                    if SegmentDict[i]["rightid"] == x:
                        rxn_line_3 = draw.Lines(
                            SegmentDict[i]["linexr"], SegmentDict[i]["lineyt"],
                            SegmentDict[i]["linexr"], SegmentDict[i]["lineyt"]+y,
                            stroke_width=line_stroke_width, # Defined in user input section RXN LINE STYLE
                            stroke = line_stroke_color, # Defined in user input section RXN LINE STYLE
                            stroke_dasharray = "10,10",
                            fill_opacity = 0, # No fill
                            close=False)
                        canvas.append(rxn_line_3)
            # Adjust line style back for future segments
            line_stroke_color=usual_line_color
            line_stroke_width=usual_line_width
            # Add reaction label to canvas if active; use global rxn label spacing to place just below bracket horizontal line
            RxnLabelHeight=SegmentDict[i]["lineyc"] - rxn_label_padding - rxn_label_font_size # Also subtract font size; text box is defined from bottom-left corner
            rxn_label = draw.Text(
                SegmentDict[i]["rxnlabel"],
                rxn_label_font_size,
                SegmentDict[i]["rxnlabelx"]+RxnLineXShift,
                RxnLabelHeight,
                fill=default_font_color)
            if show_rxn_label_placeholder == True:
                canvas.append(rxn_label)
        # Rxn line will be defined a bit differently if this is a single-segment reaction (deprotection, desulfurization, etc.)
        if "singlelinex" in SegmentDict[i]:
            rxn_line=draw.Lines(
                SegmentDict[i]["singlelinex"]+RxnLineXShift, SegmentDict[i]["singlelineyb"],
                SegmentDict[i]["singlelinex"]+RxnLineXShift, SegmentDict[i]["singlelineyt"],
                stroke_width=line_stroke_width, # Defined in user input section RXN LINE STYLE
                stroke = line_stroke_color, # Defined in user input section RXN LINE STYLE
                fill_opacity = 0, # No fill
                close=False)
            canvas.append(rxn_line)
            # Add reaction label to canvas if active; halfway between top and bottom of line
            # RxnLabelHeight=((SegmentDict[i]["singlelineyb"])+(SegmentDict[i]["singlelineyt"]))/2
            RxnLabelHeight=SegmentDict[i]["singlelineyt"] - rxn_label_padding - rxn_label_font_size # Also subtract font size; text box is defined from bottom-left corner
            rxn_label = draw.Text(
                SegmentDict[i]["rxnlabel"],
                rxn_label_font_size,
                SegmentDict[i]["rxnlabelx"]+RxnLineXShift,
                RxnLabelHeight,
                fill=default_font_color)
            if show_rxn_label_placeholder == True:
                canvas.append(rxn_label)


    # In Preview mode, return the drawSvg.drawing object; to save to file, import drawSvg module and use method canvas.saveSvg(filepath) or canvas.savePng(filepath)
    if return_canvas == True:
        return canvas

    # Save as temporary .svg and reformat tspan tags
    OutputSVGPath = f'bm_cache/temp.svg'
    if __name__ == '__main__':
        OutputSVGPath = f'output/temp.svg'
    canvas.saveSvg(OutputSVGPath)
    with open(OutputSVGPath, "r") as f:
        OutputSVGContents = f.read() # Open the file as a large string variable
    # Final reformatting step (for internal XML elements with <tag></tag>)
    OutputSVGContents = re.sub(r"\&gt;", ">", OutputSVGContents) # Replace right bracket code with >
    OutputSVGContents = re.sub(r"\&lt;", "<", OutputSVGContents) # Replace left bracket code with <
    OutputSVGContents = re.sub(r"\&amp;", "&", OutputSVGContents) # Replace ampersand code with &
    # Non-preview mode: return the formatted output file contents
    return OutputSVGContents
    ##### END MAIN FUNCTION DEF

def minusTags(InputString):
    # Returns a peptide sequence minus any ligation tags:
    SearchPattern = r"\[.*?\]"
    ReturnString = ''.join(re.split(SearchPattern, InputString))
    return(ReturnString)

def mostFrequent(InputList):
    # Determine the most common item in a list
    CountDictionary = Counter(InputList)
    return CountDictionary.most_common(1)[0][0]

def stringPercent(x):
    # Converts decimal number to string formatted as n.nn%
   return "{:.2%}".format(x)

def processStrategyQueue(internal_cys_pg='Acm',n_term_cys_pg='Thz',good_thioesters=default_good_thioesters,non_cys_thiol_penalty=2,one_pot_cys_deprotection=False,one_pot_desulfurization=False,one_pot_nterm_cys_deprotection=False,desulf_pairs={"A":"C","V":"Pen"}):
    # Primary function called repeatedly by autoBracket. Removes the first strategy from the queue list, gets all sub-strategies, and places those back in the queue list in proper order
    # To keep things straight for myself, 'strategy' is a tuple with spaces, 'bracket' is a collapsed string
    global FinalBracketList
    global StrategyQueue
    global BracketScoreDict
    # Loop through all strategies in queue, placing tags and adding results to next queue
    NextQueue=[]
    for InputTuple in StrategyQueue:
        SegmentList = InputTuple[0::2]
        LigationList = InputTuple[1::2]
        Bracket=''.join(InputTuple)
        # FIRST TIME THROUGH LOOP (i.e. first blank input strategy), this strategy will not be found in score dict. Determine which level we are placing and add to record in BracketScoreDict.
        if not Bracket in BracketScoreDict:
            BracketScoreDict.update({Bracket:{'nextlevel':getTagLevel(LigationList)}})
        Level=BracketScoreDict[Bracket]['nextlevel']
        # Loop through available ligation sites to get the range(s) in which we need to place ligation tags
        ViableRangeList = []
        PositionList = []
        for Position,Tag in enumerate(LigationList):
            # Available range = an uninterrupted stretch of blank spaces with no lower-level tags in between
            # If the space is blank, add this as a viable position
            if Tag=="[]":
                PositionList.append(Position)
            # If the tag isn't blank, and the tag contains a lower level number, this is a "split point" - mark the end of a viable range and start a new empty range
            else:
                TagNumber = re.search('(\d+)', Tag)
                if TagNumber:
                    TagNumber = TagNumber[0]
                    if int(TagNumber) < Level and len(PositionList)>0:
                        # Save position list so far (if it's not blank) to range list
                        ViableRangeList.append(PositionList)
                        # Start new blank list of viable positions
                        PositionList = []
        # At end of loop, if list of viable positions hasn't been added yet (and isn't blank), add it now
        if len(PositionList) > 0:
            ViableRangeList.append(PositionList)
        # FINISHED BRACKETS - If there are no viable positions left (i.e. no blank spaces), condense this sequence into a string and add it to our master output list of sequences
        if len(ViableRangeList) == 0:
            OutputSequence = ''.join(InputTuple)
            FinalBracketList.append(OutputSequence)
        # UNFINISHED BRACKETS - generate all combinations of tags
        else:
            # Check list of possible ranges, eliminate some options now
            UpdatedRangeList = []
            for num,ViablePositionList in enumerate(ViableRangeList):
                # If any positions will use Cys ligations, get rid of any that don't
                TempList = []
                for Position in ViablePositionList:
                    RightSegment = SegmentList[Position+1]
                    if RightSegment[0] == "C":
                        TempList.append(Position)
                if len(TempList) > 0:
                    ViablePositionList = TempList
                UpdatedRangeList.append(ViablePositionList)
            # Generate every possible permutation of viable ranges
            # Pretty nifty list comprehension trick
            PossibleCombinationList = list(itertools.product(*UpdatedRangeList))
            # PLACE TAGS AND ADD STRATEGIES TO NEXT QUEUE
            for TagCombination in PossibleCombinationList:
                StrategyIsValid=True # This will flip False if an impossible situation is detected, such as NCL of segments with an unprotected, thiolated left half (such as Pen, which can't be Thz protected)
                OutputList = list(InputTuple)
                # LOOP THROUGH TAGS IN THIS COMBINATION
                for x in TagCombination:
                    AddedSteps=1
                    CurrentTagList = [Level] # Will include explicit rxn steps as well
                    TagPosition = (x*2)+1 # Index of this tag in master segment/ligation list
                    # Current tag - evaluate if any AA substitutions need to be made
                    # First determine boundaries of this ligation; look for tags of lower levels to mark borderlines
                    # Scan middle-out; first find upper bound; by default this is the rightmost segment
                    UpperBound = len(OutputList)
                    FoundUpperBound = False
                    for i in range(TagPosition,len(OutputList),2):
                        # If we encounter a ligation tag that is a lower level than current and we haven't already set an upper bound, do so now
                        if OutputList[i] != '' and FoundUpperBound == False:
                            Tag = OutputList[i]
                            TagNumber = re.search('(\d+)', Tag)
                            if TagNumber:
                                TagNumber = TagNumber[0]
                                if int(TagNumber) <= Level:
                                    UpperBound = i # Index just past rightmost segment
                                    FoundUpperBound = True
                    # Do the same to find lower bound; by default this is the leftmost segment
                    LowerBound = 0
                    FoundLowerBound = False
                    # Go in reverse order, from tag position to 0 backwards
                    for i in sorted(range(1,TagPosition,2),reverse=True):
                        if OutputList[i] != '' and FoundLowerBound == False:
                            Tag = OutputList[i]
                            TagNumber = re.search('(\d+)', Tag)
                            if TagNumber:
                                TagNumber = TagNumber[0]
                                if int(TagNumber) < Level:
                                    LowerBound = i+1 # Index of leftmost segment
                                    FoundLowerBound = True
                    # MAKE SUBSTITUTIONS
                    # If this is a non-Cys ligation, it means that ALL empty spaces within these boundaries are non-Cys ligations
                    RightSegment = OutputList[TagPosition+1]
                    # These toggles indicate whether we need to include the keyword or not, after making substitutions
                    NeedsDesulfStep=False
                    CysWasProtected=False
                    NTermCysWasProtected=False
                    if RightSegment[0] in desulf_pairs: # Pairs are defined in input variable, a dictionary of After:Before
                        SegmentPositionsList=list(range(LowerBound,UpperBound,2)) # Indices of all relevant segments; UpperBound should be the index just past rightmost segment
                        # First, protect any internal Cys in these segments, if present; mark if I had to do that
                        for SegmentPosition in SegmentPositionsList:
                            Segment=OutputList[SegmentPosition]
                            NewSegment=subAA('C',f'C-{internal_cys_pg}',Segment)
                            if Segment!=NewSegment:
                                # Update segment in list, and record that we have done a Cys protection here
                                OutputList[SegmentPosition]=NewSegment
                                CysWasProtected=True
                        # If any internal Cys were protected, and there is an N-terminal Cys on the N-terminal segment, protect that the same way here (rather than using unique N-term Cys protecting group in a separate step)
                        if CysWasProtected==True:
                            LeftmostSegment=OutputList[LowerBound]
                            LeftmostSegmentSubbed=subAA('C',f'C-{internal_cys_pg}',LeftmostSegment,FirstOnly=True)
                            if LeftmostSegmentSubbed!=LeftmostSegment:
                                OutputList[LowerBound]=LeftmostSegmentSubbed
                        # Then, substitute all N-terminal AA with the proper desulfurization pair (C>A), except the leftmost segment
                        for SegmentPosition in SegmentPositionsList:
                            if SegmentPosition!=LowerBound:
                                Segment=OutputList[SegmentPosition]
                                # Substitute first AA for proper keyword depending on AA identity
                                NewSegment=Segment
                                for after_key in desulf_pairs:
                                    before_key=desulf_pairs[after_key]
                                    SubString=f'{before_key}>{after_key}'
                                    NewSegment=subAA(after_key,SubString,NewSegment,FirstOnly=True)
                                if NewSegment!=Segment:
                                    OutputList[SegmentPosition]=NewSegment
                                    NeedsDesulfStep=True
                        #/Ala or Val substitutions
                    # Leftmost segment only - protect N-terminal Cys if it exists and hasn't already been protected from Ala desulf
                    if CysWasProtected==False:
                        LeftmostSegment=OutputList[LowerBound]
                        LeftmostSegmentSubbed=subAA('C',f'C-{n_term_cys_pg}',LeftmostSegment,FirstOnly=True)
                        if LeftmostSegmentSubbed!=LeftmostSegment:
                            OutputList[LowerBound]=LeftmostSegmentSubbed
                            NTermCysWasProtected=True
                    # CHECK FOR SELF-CYCLIZATION; only with expanded thiols beyond Cys and Ala, which are considered incompatible with Thz or similar protection methods
                    LeftmostSegment=OutputList[LowerBound]
                    if LeftmostSegment[0]=="(":
                        LeftmostAA=re.split(r'[\(\)>]',LeftmostSegment)[1]
                        if LeftmostAA in desulf_pairs.values() and LeftmostAA!="C":
                            StrategyIsValid=False
                    # UPDATE TAG KEYWORDS
                    if NeedsDesulfStep==True:
                        if one_pot_desulfurization==True:
                            CurrentTagList.append('Desulf')
                        else:
                            CurrentTagList.append('Desulf*')
                    if CysWasProtected==True:
                        if one_pot_cys_deprotection==True:
                            CurrentTagList.append(internal_cys_pg)
                        else:
                            CurrentTagList.append(internal_cys_pg+'*')
                    if NTermCysWasProtected==True:
                        if one_pot_nterm_cys_deprotection==True:
                            CurrentTagList.append(n_term_cys_pg)
                        else:
                            CurrentTagList.append(n_term_cys_pg+'*')
                    # FORMAT TAG AS STRING AND ADD IN PLACE
                    CurrentTagList = [str(x) for x in CurrentTagList]
                    CurrentTag='['
                    for n,TagText in enumerate(CurrentTagList):
                        if TagText[-1]=='*': # Asterisk indicates explicit step in label; change to [1][step] syntax
                            CurrentTag+=(']['+TagText.rstrip('*'))
                        else:
                            if n>0:
                                CurrentTag+=','
                            CurrentTag+=TagText
                    CurrentTag+=']'
                    OutputList[TagPosition] = CurrentTag
                    #/loop through tags in tag combination
                # SCORE STRATEGY AND ADD TO NEXT QUEUE
                if StrategyIsValid==True:
                    OutputSequence = ''.join(OutputList)
                    BracketScoreDict.update({OutputSequence:scoreBracket(OutputSequence,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty)})
                    # Add to next queue, and record the next level to check for this strategy
                    OutputTuple = tuple(OutputList)
                    NextQueue.append(OutputTuple)
                    BracketScoreDict[OutputSequence].update({'nextlevel':Level+1})
                    #/loop through tag combinations
            #/get tag combinations for unfinished brackets
        #/loop through strategies in queue
    # CLEAR STRATEGY QUEUE, SORT AND UPDATE NEXT QUEUE
    SortedNextQueue = sorted(NextQueue, key = lambda x: countTags(x),reverse=True) # Highest number of tags placed, to lowest
    # If there are any complete strategies in the queue, add them to the final strategy list now (faster to do this at end of loop then as part of next loop)
    StrategiesToRemove=[]
    for Strategy in SortedNextQueue:
        if not '[]' in Strategy:
            FinalBracketList.append(''.join(Strategy))
            StrategiesToRemove.append(Strategy)
    for Strategy in StrategiesToRemove:
        SortedNextQueue.remove(Strategy)
    StrategyQueue=SortedNextQueue[:]


def recursiveRange(InputString, CurrentLevel, i=0):
    # Part of main function; takes input string (XXXXXXX + ligation tags) and converts to nested list format with range of each segment
    # i is the index of the first AA (in sequence list)
    # This function converts the input string to a nested list form that is easier for Python to interpret, not necessarily the user.
    # Recursion will stop if it doesn't detect the current level in the sequence.
    # Break sequence up into a list of numbers:
    SplitAtNumbers = re.split(r"(\d+)",InputString)
    CurrentLevelCount = SplitAtNumbers.count(str(CurrentLevel))
    if CurrentLevelCount == 0:
        # ERROR CHECK: If current level isn't found but other ligation tags are, display the following error.
        LigationTagSearchPattern = r"\[.*?\]"
        ListOfLigations = re.findall(LigationTagSearchPattern,InputString)
        NumberOfLigations = len(ListOfLigations)
        if NumberOfLigations > 0:
            print(f"Error: Encountered ligation tag(s) {ListOfLigations} in InputString {InputString}, but no tags with CurrentLevel = {CurrentLevel}. Ensure no levels were skipped.")
            sys.exit(f"Error: Encountered ligation tag(s) {ListOfLigations} in InputString {InputString}, but no tags with CurrentLevel = {CurrentLevel}. Ensure no levels were skipped.")
        # If no ligation tags are found, this a terminal segment and the current iteration of the fxn ends. Return a range of the boundaries of the current segment (based on index of sequence list, not AA number)
        else:
            return range(i,i+len(InputString))
    else:
        # If I formulate my search pattern with capturing (parentheses), it should return a list with exactly 3 items (if this is all done correctly)
        SearchPattern = r"(\["+str(CurrentLevel)+r".*?\])"
        list_to_return = re.split(SearchPattern, InputString)
        # Second item in list (ligation tag) needs to be reformatted to a list, not a string.
        list_to_return[1] = stringToList(list_to_return[1])
        # Error checking: If current level is found twice, display the following error.
        if len(list_to_return)>3:
            print(f'Error: Encountered two adjacent [{CurrentLevel}] tags; must be separated by a tag of [{CurrentLevel-1}] or lower.')
            sys.exit(f'Error: Encountered two adjacent [{CurrentLevel}] tags; must be separated by a tag of [{CurrentLevel-1}] or lower.')
        # There should be exactly 3 items in the list if done correctly. Let's give them names.
        LeftHalfString = list_to_return[0]
        LigationTag = list_to_return[1]
        RightHalfString = list_to_return[2]
        # // BEGIN RECURSION
        # Update the segment list to convert L and R halves to their own sub-lists of segments and ligations. Start right half at first AA of right segment
        list_to_return[0] = recursiveRange(LeftHalfString, CurrentLevel + 1, i)
        list_to_return[2] = recursiveRange(RightHalfString, CurrentLevel + 1, i+len(''.join(re.split(r"\[.*?\]",LeftHalfString))))
        # // END RECURSION
        return list_to_return

def relativeCanvasPosition(AAIndex):
    # Converts between "position in protein" and "position on canvas" relative to bottom segment width
    return AAIndex / len(SequenceList) * bottom_segment_width


def scoreBracket(InputSequence,good_thioesters=default_good_thioesters,non_cys_thiol_penalty=2,typical_rxn_yield=0.5):
    # "Top down" scoring function; input is a partial or complete bracket string
    # Remove any empty tags [], replacing with special character $ (easier for path lengths fxn to handle)
    InputSequence="$".join(InputSequence.split("[]"))
    # Returns a dictionary with various scoring metrics
    Score={}
    # SINGLE SEGMENTS: return a blank score dictionary
    if not '[' in InputSequence:
        Score.update({
        'maxpath':0,
        'sumpath':0,
        'avgpath':0,
        'thioester':0,
        'minyield':0,
        'avgyield':0,
        'steps':0,
        'segments':0
        })
        return Score
    # Split into list of segments & ligations
    LigationTagSearchPattern = r"(\[.*?\])"
    SegmentLigationList = re.split(LigationTagSearchPattern,InputSequence)
    # If there is a blank segment, this is actually a blank space between two explicit rxns; reformat as a keyword with asterisk
    while '' in SegmentLigationList:
        i=SegmentLigationList.index('')
        FirstTag=stringToList(SegmentLigationList[i-1])
        SecondTag=stringToList(SegmentLigationList[i+1])
        ExplicitStepLabel = SecondTag[0]+'*'
        SecondTag[0]=ExplicitStepLabel
        CombinedTag=listToString(FirstTag+SecondTag)
        TempList = SegmentLigationList[0:i-1]+[CombinedTag]+SegmentLigationList[i+2:]
        SegmentLigationList = TempList[:]
    InputSequence=''.join(SegmentLigationList)
    # Determine path length of each segment and add scores
    AllPathLengths = getPathLengths(InputSequence,good_thioesters=good_thioesters,non_cys_thiol_penalty=non_cys_thiol_penalty)
    Score.update({
        'maxpath':max(AllPathLengths.values()),
        'sumpath':sum(AllPathLengths.values()),
        'avgpath':sum(AllPathLengths.values())/len(AllPathLengths),
        'segments':len(AllPathLengths)
        })
    # # Determine contribution from thioester penalty...note, disabled by default because it doubles the calculation time of scoring brackets
    # IdealPathLengths=getPathLengths(InputSequence,good_thioesters=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ)'),non_cys_thiol_penalty=non_cys_thiol_penalty)
    # ThioesterPenalty=sum(AllPathLengths.values())-sum(IdealPathLengths.values())
    # Score.update({
    #     'thioester':ThioesterPenalty
    #     })
    # Temporarily put this in place so that Auto-Bracket still works
    Score.update({
        'thioester':0
        })
    # Convert to segment yields and add scores
    AllYields={}
    for Segment in AllPathLengths:
        PathLength=AllPathLengths[Segment]
        Yield=typical_rxn_yield**PathLength
        AllYields.update({Segment:Yield})
    Score.update({
        'minyield':min(AllYields.values()),
        'avgyield':sum(AllYields.values())/len(AllYields)
        })
    # Add score for total reaction steps; only count ligations & explicit steps
    LigationList=SegmentLigationList[1::2]
    TotalSteps=0
    for Tag in LigationList:
        TotalSteps+=1
        TotalSteps+=Tag.count('*')
    Score.update({
        'steps':TotalSteps
        })
    return Score


def split(InputString):
    # Splits a text string into a list of characters
    return [char for char in InputString]


def stringToList(InputString):
    # Converts a string that looks like a list (e.g. "[NCL,0]") to a correctly-interpreted list object.
    NewList = InputString.strip('][').split(',')
    return NewList


def subAA(SingleLetter,Substitution,InputString,FirstOnly=False,DoneSubbing=False):
        # Takes input in the form of a single segment, and outputs the segment with all AA, or only the first AA, substituted with the specified string
        # Substitutions must be exactly what you want replaced in the string, not including parentheses
        InsideParentheses = False
        CharList = []
        for Char in InputString:
            if DoneSubbing == True: # If requested to only sub the first AA and we've already done that, this will be True
                CharList.append(Char)
            # If we're not done subbing yet, evaluate whether to make the substitution here
            else:
                # If we encounter an opening parentheses, begin capturing all contents in a separate list
                if Char == "(":
                    InsideParentheses = True
                    InsideParenthesesList = []
                if InsideParentheses == True:
                    # Looking for pattern (C>A) to replace with (C-Acm>A)
                    # If character is > and only other character seen so far is an unprotected SingleLetter, replace it with the substitution
                    if Char == '>' and len(InsideParenthesesList) == 2 and InsideParenthesesList[1] == SingleLetter:
                        InsideParenthesesList[1] = Substitution
                    # Regardless, add all characters to InsideParenthesesList
                    InsideParenthesesList.append(Char)
                    # When a closing parentheses is encountered, condense into a single string
                    if Char == ')':
                        InsideParentheses = False
                        Char = ''.join(InsideParenthesesList)
                if InsideParentheses == False and Char == SingleLetter:
                    # Change this character to its substitution and put parentheses around it
                    Char = f'({Substitution})'
                if InsideParentheses == False:
                    # Add character (or short phrase) to ongoing string
                    CharList.append(Char)
                    # If it's requested to only sub the first character, toggle now
                    if FirstOnly == True:
                        DoneSubbing = True
        # After looping through all characters, return the final condensed string
        return ''.join(CharList)

def subNumbers(InputString,runfxn=True,FontSize=0):
    # Convert all numbers in a text block to subscript (to be passed to SVG file within <text></text> tag)
    # Works specifically with drawSvg module and drawBracket fxn
    # Ignore for PNG, returning original string
    if runfxn==False:
        return InputString
    if FontSize==0:
        FontSize=segment_label_font_size
    # Split input string by numbers; include commas as well
    SplitAtNumbers = re.split(r"([\d,]+)", InputString)
    # Define start and end tags for svg XML code
    TspanStart = r'<tspan baseline-shift = "sub" style="font-size: ' + str(round(FontSize*10/12,0)) + r';">'
    TspanEnd = r'</tspan>'
    # Make an alternating list of appropriate length
    Length = int(len(SplitAtNumbers)/2)
    TspanList = [TspanStart,TspanEnd]*Length
    # Combine the two lists in alternating fashion
    CombinedList = [None]*(len(SplitAtNumbers)+len(TspanList))
    CombinedList[::2] = SplitAtNumbers
    CombinedList[1::2] = TspanList
    # Converge to string and return
    return ''.join(CombinedList)



#-------------------------------------------------------------------
#--------------RUNTIME START----------------------------------------
#-------------------------------------------------------------------
if __name__ == '__main__':
    gettime('Defined all functions')
    timestamp=datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Below are blocks which may be un-commented to test changes in the main functions

    # ##### BRACKET DRAWING TEST
    # # Update the InputSequence with custom bracket notation
    # drawingtestsequence ='(Nle)IEKLRNIAIIAHVDHGKTTLVDKLLQQSGTFDSR[4](C>A)ETQERV(Nle)DSNDLEKERGITILAKNTAIKWNDYRINIVDTPGH[3](C-TfaThz>A)DFGGEVERV(Nle)S(Nle)VDSVLLVVDAFDGP(Nle)PQTRFVTKKAF[4,TfaThz](C>A)YGLKPIVVINKVDRPG[2][Desulf][Acm](C-TfaThz>A)RPDWVVDQVFDLFVNLDATDEQLDFPIVYASALNGIAGLDHED(Nle)[4](C>A)EDMTPLYQAIVDHVPAPDVDLDGPFQMQISQLDYNSYVG[3,TfaThz](Pen>V)IGIGRIKRGKVKPNQQVTIIDSEGKTRNAKVGKVLGHLGLERIETDL[4](C>A)EAGDIVAITGLGELNISDTV(C-Acm)DTQNVEALPALSVDEPTVSMFF[1](C-Acm)VNTSPF(C-Acm)GKEGKFVTSRQILDRLNKELVHNV[3](C-TfaThz>A)LRVEETEDADAFRVSGRGELHLSVLIENMRREGFELA[4,TfaThz](Pen>V)SRPKVIFREIDGRKQEPYENVTLDVEEQHQGSVMQALGERKGDLKNMNPDGKGRVRLDY[2][Desulf][Acm](Pen>V)IPSRGLIGFRSEFMTMTSGTGLLYSTFSHYDDVRPGEVGQRQNG[4](Pen>V)LISNGQGKAVAFALFGLQDRGKLFLGHG[3](C-TfaThz>A)EVYEGQIIGIHSRSNDLTVN(C-Acm)LTGKKLTNMR[4,TfaThz](C>A)SGTDEAVVLVPPIRMTLEQALEFIDDDELVEVTPTSIRIRKRHLTENDRRRANRAPKDD' # Example of E. coli protein BipA; complex bracket with basically every feature, good for testing changes
    # with open(f'output/test_{timestamp}.svg', 'w') as f:
    #     f.write(makeBracketFigure(drawingtestsequence)) # Change fxn variables here; otherwise defaults will be used
    # gettime('BracketMaker Test Complete')


    # ##### BRACKET SCORING TEST
    # # Update with bracket notation input
    # scoretestsequence ='(Nle)IEKLRNIAIIAHVDHGKTTLVDKLLQQSGTFDSR[4](C>A)ETQERV(Nle)DSNDLEKERGITILAKNTAIKWNDYRINIVDTPGH[3](C-TfaThz>A)DFGGEVERV(Nle)S(Nle)VDSVLLVVDAFDGP(Nle)PQTRFVTKKAF[4,TfaThz](C>A)YGLKPIVVINKVDRPG[2][Desulf][Acm](C-TfaThz>A)RPDWVVDQVFDLFVNLDATDEQLDFPIVYASALNGIAGLDHED(Nle)[4](C>A)EDMTPLYQAIVDHVPAPDVDLDGPFQMQISQLDYNSYVG[3,TfaThz](Pen>V)IGIGRIKRGKVKPNQQVTIIDSEGKTRNAKVGKVLGHLGLERIETDL[4](C>A)EAGDIVAITGLGELNISDTV(C-Acm)DTQNVEALPALSVDEPTVSMFF[1](C-Acm)VNTSPF(C-Acm)GKEGKFVTSRQILDRLNKELVHNV[3](C-TfaThz>A)LRVEETEDADAFRVSGRGELHLSVLIENMRREGFELA[4,TfaThz](Pen>V)SRPKVIFREIDGRKQEPYENVTLDVEEQHQGSVMQALGERKGDLKNMNPDGKGRVRLDY[2][Desulf][Acm](Pen>V)IPSRGLIGFRSEFMTMTSGTGLLYSTFSHYDDVRPGEVGQRQNG[4](Pen>V)LISNGQGKAVAFALFGLQDRGKLFLGHG[3](C-TfaThz>A)EVYEGQIIGIHSRSNDLTVN(C-Acm)LTGKKLTNMR[4,TfaThz](C>A)SGTDEAVVLVPPIRMTLEQALEFIDDDELVEVTPTSIRIRKRHLTENDRRRANRAPKDD' # Example of E. coli protein BipA; complex bracket with basically every feature, good for testing changes
    # print(f'Test Bracket Scores:')
    # print(scoreBracket(scoretestsequence))


    # ##### AUTOBRACKET TEST
    # # Update with empty bracket - place [] between segments
    # autobrackettestsequence="MIEKLRNIAIIAHVDHGKTTLVDKLLQQSGTFDSR[]AETQERVMDSNDLEKERGITILAKNTAIKWNDYRINIVDTPGH[]ADFGGEVERVMSMVDSVLLVVDAFDGPMPQTRFVTKKAF[]AYGLKPIVVINKVDRPG[]ARPDWVVDQVFDLFVNLDATDEQLDFPIVYASALNGIAGLDHEDM[]AEDMTPLYQAIVDHVPAPDVDLDGPFQMQISQLDYNSYVG[]VIGIGRIKRGKVKPNQQVTIIDSEGKTRNAKVGKVLGHLGLERIETDL[]AEAGDIVAITGLGELNISDTVCDTQNVEALPALSVDEPTVSMFF[]CVNTSPFCGKEGKFVTSRQILDRLNKELVHNV[]ALRVEETEDADAFRVSGRGELHLSVLIENMRREGFELA[]VSRPKVIFREIDGRKQEPYENVTLDVEEQHQGSVMQALGERKGDLKNMNPDGKGRVRLDY[]VIPSRGLIGFRSEFMTMTSGTGLLYSTFSHYDDVRPGEVGQRQNG[]VLISNGQGKAVAFALFGLQDRGKLFLGHG[]AEVYEGQIIGIHSRSNDLTVNCLTGKKLTNMR[]ASGTDEAVVLVPPIRMTLEQALEFIDDDELVEVTPTSIRIRKRHLTENDRRRANRAPKDD" # 50S ribosomal subunit assembly factor BipA; requires Val ligations; winning strategy has a mix of Val, Ala and Cys
    # autoBracket(autobrackettestsequence,output_scores_file=True,filename=f'AutoBracket_{timestamp}',returnnumber=100,verbose=True) # Add other parameters if desired
    # gettime('AutoBracket Test Complete')
