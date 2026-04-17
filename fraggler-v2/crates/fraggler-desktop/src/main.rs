use std::path::PathBuf;

use anyhow::Result;
use slint::{ComponentHandle, SharedString};

slint::include_modules!();

fn main() -> Result<()> {
    let app = MainWindow::new()?;
    let readme_path = workspace_readme_path();

    {
        let weak = app.as_weak();
        app.on_show_contract(move || {
            if let Some(app) = weak.upgrade() {
                let log = format!(
                    "{}\n[contract]\n- v1 JSON contract scaffolded\n- CLI emits JSON lines\n- Rust core currently returns NotImplemented until engine port begins\n",
                    app.get_log_text()
                );
                app.set_log_text(SharedString::from(log));
                app.set_status_text(SharedString::from(
                    "Contract scaffold is ready. Next step is engine porting.",
                ));
            }
        });
    }

    {
        let weak = app.as_weak();
        app.on_run_demo(move || {
            if let Some(app) = weak.upgrade() {
                let log = format!(
                    "{}\n[desktop]\nDesktop shell scaffold is live, but it is not yet connected to a real Rust engine.\n",
                    app.get_log_text()
                );
                app.set_log_text(SharedString::from(log));
                app.set_status_text(SharedString::from(
                    "Desktop shell scaffolded. Engine integration is still pending.",
                ));
            }
        });
    }

    app.on_open_docs(move || {
        if let Err(err) = open::that(&readme_path) {
            eprintln!("failed to open docs: {err}");
        }
    });

    app.run()?;
    Ok(())
}

fn workspace_readme_path() -> PathBuf {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .map(|path| path.join("README.md"))
        .unwrap_or_else(|| PathBuf::from("README.md"))
}
