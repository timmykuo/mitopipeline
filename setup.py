from setuptools import setup

setup(
    name="mitopipeline",
    version="0.13",
    packages=['mitopipeline'],
    py_modules=['mitopipeline'],

    #add samtools, bwa, etc. here
    install_requires=[
        'jinja2',
        'luigi',
        'setuptools',
    ],

    setup_requires=[
        'jinja2',
        'luigi',
        'setuptools'
    ],

    entry_points = '''
        [console_scripts]
        mitopipeline=mitopipeline.cmdline:run
    ''',
    package_data={'mitopipeline': ['/mitopipeline/steps/*']},

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
