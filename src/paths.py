from pathlib import Path

def find_project_root(start: Path) -> Path:
    """Sube carpetas hasta encontrar .git o README.md."""
    start = start.resolve()
    for p in [start] + list(start.parents):
        if (p / ".git").exists() or (p / "README.md").exists():
            return p
    raise RuntimeError(
        "No se pudo detectar PROJECT_ROOT. Ejecuta el notebook desde dentro del repo "
        "(o asegúrate de que existe README.md en la raíz)."
    )

def project_paths(cwd: Path) -> dict:
    root = find_project_root(cwd)
    return {
        "NOTEBOOK_DIR": cwd.resolve(),
        "PROJECT_ROOT": root,
        "CONFIG_DIR": root / "config",
        "DATA_DIR": root / "data",
        "RESULTS_DIR": root / "results",
        "FIGURES_DIR": root / "figures",
        "SRC_DIR": root / "src",
        "R_DIR": root / "R",
        "NOTEBOOKS_DIR": root / "notebooks",
    }
