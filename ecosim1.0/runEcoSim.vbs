Dim args
Dim command
Dim debug
Dim path
Dim shell

' Set the command for EcoSim.
command = "java -jar bin\EcoSim.jar"

' Grab the path.
path = Replace (Wscript.ScriptFullName, Wscript.ScriptName, "")

' Create a shell object.
Set shell = Wscript.CreateObject ("Wscript.Shell")
	
' Change to the current working directory.
shell.CurrentDirectory = path

' Check the arguments to see if the debug switch is on.
debug = 0
If Wscript.Arguments.Count > 0 Then
  For I = 0 to WScript.Arguments.Count - 1
    args = args & " " & Wscript.Arguments.Item(I)
  Next
  For I = 0 to WScript.Arguments.Count - 1
    If Wscript.Arguments.Item(I) = "--debug" Then
      debug = 1  
    ElseIf Wscript.Arguments.Item(I) = "-d" Then
      debug = 1
    End If
  Next
End If

' Run EcoSim.
shell.Run command & " " & args, debug, FALSE

