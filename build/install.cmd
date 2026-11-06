@echo off
setlocal
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1"
set "install_status=%ERRORLEVEL%"
if not "%install_status%"=="0" echo Installation failed. See the error above.
pause
exit /b %install_status%
