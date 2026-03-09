import PyInstaller.__main__
import sys
import os

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
