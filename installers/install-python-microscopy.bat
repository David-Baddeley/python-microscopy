@echo off
setlocal EnableDelayedExpansion

:: ---- Config (mirrors installer_defines.sh — keep in sync) ----
set "TARGET_PYTHON=3.13"
set "PACKAGE_NAME=python-microscopy"
set "ENTRY_POINTS=PYMEAcquire PYMEImage PYMEVis PYMEClusterOfOne"
set "DEFAULT_DEST=%USERPROFILE%\PYME"

:: ---- Destination (first positional arg, else default) ----
if not "%~1"=="" (set "DEST=%~1") else (set "DEST=%DEFAULT_DEST%")
echo Installing PYME to: !DEST!
if not exist "!DEST!\" mkdir "!DEST!"

:: ---- Locate or download uv ----
:: Context A (CI): uv is in PATH via astral-sh/setup-uv action — the where check passes immediately.
:: Context B (end-user): uv.exe is downloaded to DEST\bin\ and used from there.
where uv >nul 2>&1
if not errorlevel 1 (
    set "UV=uv"
) else (
    call :download_uv
    if errorlevel 1 exit /b 1
)

:: ---- Standalone Python installation via py-app-standalone ----
echo Creating standalone Python installation...
"!UV!" tool run py-app-standalone --python-version !TARGET_PYTHON! --target "!DEST!" !PACKAGE_NAME!
if errorlevel 1 (echo ERROR: py-app-standalone failed & exit /b 1)

:: ---- Normalize cpython-* to a stable directory name ----
for /d %%d in ("!DEST!\cpython-*") do move "%%d" "!DEST!\python" >nul
if not exist "!DEST!\python\" (echo ERROR: standalone Python directory not found & exit /b 1)

:: ---- Entry point .cmd wrappers ----
for %%e in (%ENTRY_POINTS%) do call :mk_wrapper "%%e"

:: ---- Console helper (adds python dir to PATH for the session) ----
(
    echo @echo off
    echo set "PATH=%%~dp0python;%%~dp0python\Scripts;%%PATH%%"
    echo cmd.exe /k
) > "!DEST!\pyme-console.cmd"

echo.
echo Installation complete.
echo   Wrappers: %ENTRY_POINTS%
echo   Add !DEST! to your PATH, or run the .cmd files directly.
echo   Activated console: !DEST!\pyme-console.cmd
goto :eof


:: ----------------------------------------------------------------
:download_uv
:: Downloads uv.exe (x86-64) into DEST\bin\ via curl + PowerShell.
:: Requires Windows 10 1803+ (curl.exe and tar.exe built-in).
if not exist "!DEST!\bin\" mkdir "!DEST!\bin"
set "UV_ZIP=%TEMP%\uv_download.zip"
echo Downloading uv (x86-64)...
curl -fsSLo "%UV_ZIP%" "https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-pc-windows-msvc.zip"
if errorlevel 1 (echo ERROR: Failed to download uv & exit /b 1)
powershell -NoProfile -Command "Expand-Archive -LiteralPath '%UV_ZIP%' -DestinationPath '!DEST!\bin' -Force"
if errorlevel 1 (echo ERROR: Failed to extract uv & exit /b 1)
del "%UV_ZIP%"
set "UV=!DEST!\bin\uv.exe"
exit /b 0


:: ----------------------------------------------------------------
:mk_wrapper
:: Uses %~dp0 so wrappers stay valid if the install folder is moved.
set "_EP=%~1"
(
    echo @echo off
    echo "%%~dp0python\Scripts\!_EP!.exe" %%*
) > "!DEST!\!_EP!.cmd"
exit /b 0
