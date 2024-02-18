DIM objShell
set objShell=wscript.createObject("wscript.shell")
iReturn=objShell.Run("dbgrn.exe /D:\dbgrn.exe", 0, TRUE)
