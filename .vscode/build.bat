@echo off
REM Usage: build.bat <source.cpp> <output.exe> [run]
REM Loads the MSVC environment, then compiles the given source file.
setlocal

REM Locate the VC toolchain with vswhere so this survives Visual Studio updates
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSINSTALL="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSINSTALL=%%i"
)
if not defined VSINSTALL set "VSINSTALL=C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools"

set "VCVARS=%VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VCVARS%" (
  echo ERROR: vcvars64.bat not found at "%VCVARS%"
  exit /b 1
)

REM vcvars prints a harmless 'vswhere.exe not recognized' line on Build Tools installs
call "%VCVARS%" >nul 2>&1

set "OBJDIR=%~dp0..\build"
if not exist "%OBJDIR%" mkdir "%OBJDIR%"

REM Eigen is header-only; point at the extracted release (see README / vcxproj)
set "EIGEN_DIR=%USERPROFILE%\eigen-3.4.0"
if not exist "%EIGEN_DIR%\Eigen\Dense" (
  echo ERROR: Eigen not found at "%EIGEN_DIR%"
  echo Download https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip and extract it there.
  exit /b 1
)

REM _USE_MATH_DEFINES is set here, not in a header: it only has an effect before
REM the first <cmath> in a translation unit, so a header-level #define silently
REM stops working as soon as the include order changes.
cl /nologo /EHsc /std:c++17 /O2 /DNDEBUG /Zi /W3 ^
   /D_USE_MATH_DEFINES ^
   /I "%EIGEN_DIR%" ^
   /Fe:"%~2" /Fo:"%OBJDIR%\\" /Fd:"%OBJDIR%\\" ^
   "%~1"
if errorlevel 1 exit /b 1

if /i "%~3"=="run" (
  echo.
  "%~2"
)

exit /b 0
