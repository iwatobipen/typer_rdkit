from setuptools import setup

setup(
    name="rdkit-cli",
    version="0.1",
    py_modules=[],
    install_requires=[
        "typer",
        "typer-cli",
        "useful_rdkit_utils",
        "numpy",
        "rdkit"
    ],
    entry_points="""
    [console_scripts]
    rdkit_cli=rdkit_cli:main
    """,
)
