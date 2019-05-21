from distutils.core import setup

setup(
    name="scrap",
    version="0.1",
    packages=[
        "scrap",
        "scrap.tags",
        "scrap.chromium"
    ],
    install_requires=[
        "pandas",
        "pepars"
    ]
)
