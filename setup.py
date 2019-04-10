from distutils.core import setup

setup(
    name="mitopipeline",
    version="0.1",
    packages=['mitopipeline'],

    #add samtools, bwa, etc. here
    install_requires=[
        'jinja2'
    ],

    package_data={},

    #metadata
    author="Timothy Kuo",
    author_email="timykuo@gmail.com",
    description="Mitochondrial Genome Pipeline",
    license="MIT",
    keywords="genetics research data pipeline dependency management mitochondrial genome",
    url="https://github.com/timmykuo/mitopipeline",
    download_url="",

    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    project_urls={
        "Bug Tracker": "",
        "Documentation": "",
        "Source Code": "",
    }
)
