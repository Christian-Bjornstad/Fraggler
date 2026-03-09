import PyInstaller.__main__
import sys
import os
import dis

# Monkey-patch dis._get_const_info to swallow known Python 3.10 bytecode IndexErrors
_orig_get_const_info = getattr(dis, '_get_const_info', None)
if _orig_get_const_info:
    def _patched_get_const_info(arg, constants):
        try:
            return _orig_get_const_info(arg, constants)
        except IndexError:
            return arg, repr(arg)
    dis._get_const_info = _patched_get_const_info

def build_app():
    print(f"Building Fraggler Diagnostics for {sys.platform}...")
    
    args = [
        'qt_app.py',
        '--name=Fraggler',
        '--noconfirm',
        '--windowed', # No console window
        '--clean',
        '--add-data=assets:assets', # Include CSS / images if needed
        '--hidden-import=PyQt6',
        '--hidden-import=pandas',
        '--hidden-import=plotly',
    ]
    
    if sys.platform == 'darwin':
        args.append('--icon=assets/app_icon.icns')
        # macOS specific bundle settings can be added here
    elif sys.platform == 'win32':
        args.append('--icon=assets/app_icon.ico')
    elif sys.platform == 'linux':
        pass
        
    PyInstaller.__main__.run(args)
    print("\nBuild complete! Check the 'dist' directory.")

if __name__ == "__main__":
    build_app()
