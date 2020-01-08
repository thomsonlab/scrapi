from distutils.core import setup

setup(
    name="scrap",
    version="0.3.1",
    packages=[
        "scrap",
        "scrap.tags",
        "scrap.chromium",
        "scrap.plotting",
        "scrap.dataset",
        "scrap.utils"
    ],
    entry_points={
        "console_scripts": [
            "scrap-preprocess=scrap.dataset.preprocess:preprocess",
            "scrap-init-ws=scrap.dataset.initialize:initialize"
        ]
    },
    install_requires=[
        "pandas",
        "pepars @ git+https://github.com/gradinarulab/sparsedat-py.git",
        "sparsedat @ git+https://github.com/thomsonlab/sparsedat-py.git",
        "numpy",
        "plotly",
        "h5py",
        "sklearn",
        "scipy",
        "statsmodels"
    ]
)
