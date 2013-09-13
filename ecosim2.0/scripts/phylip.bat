@echo off
:: @param 1 The Phylip program to run.
:: @param 2 The directory containing the infile.
:: @param 3 The debug option.

:: Change to the directory containing the infile.
cd "%2"

:: Remove the outfile and outtree files if they exist.
if exist outfile del outfile
if exist outtree del outtree

:: Look for the debug option.
if "%3"=="true" (
  :: Run the Phylip program with debugging on.
  type input | "%1"
) else (
  :: Run the Phylip program with debugging off.
  type input | "%1" > screenout
  del screenout
)
