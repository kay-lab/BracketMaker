# BracketMaker
BracketMaker, a graphical workspace for chemical protein synthesis - make bracket figures and predict synthetic routes.

V1 beta for early access, feedback welcome! Official release will accompany publication of manuscript (in review).

Executables for MacOS and Windows may be found in the Releases tab. See User Manual PDF for instructions on installation and use.
Please note that the Windows version was prepared on Windows 10 and may only work on Windows 10 and above. The MacOS version was prepared on MacOS 10.13 (High Sierra) and may only work with OS X versions 10.13 and above. The MacOS executable will only work if the the *cairo* graphics library is installed: https://cairographics.org/ (see User Manual for instructions on installing all prerequisites).

----

## Program Description

The bracketmaker_gui.py script encodes all elements of the tkinter-based GUI. The novel functions of BracketMaker used by this GUI are defined in scripts/bracketmaker.py, including text-to-figure image generation, interpretation and scoring of brackets, and automatic bracket assembly from a list of segments. See our BracketMaker publication (pending) for more information on these functions and the development of the program.

To generate a bracket figure, type a protein sequence into the main text box of the GUI, and insert numbered ligation junctions [#] between segments. See the User Manual for a full description of syntax and a bracket-building tutorial.

BracketMaker also includes an AutoBracket feature ideal for evaluating new protein targets. Aligator (see our /kay-lab/Aligator repository) and AutoBracket may be used in combination to quickly generate a list of synthetic strategies for an input protein sequence. See User Manual for more information.

This script is compatible with Python 3.10 and higher. If the user is working with the source code, *tkinter*, *importlib*, and *drawSvg* packages must be installed, as well as the *cairo* vector graphics library: https://cairographics.org/ .
