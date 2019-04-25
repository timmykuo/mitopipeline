from setuptools import setup

setup(
    name="mitopipeline",
    version="1.0",
    packages=['mitopipeline'],
    package_dir={'mitopipeline': 'mitopipeline'},
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
    package_data={'mitopipeline': ['steps/*', 'steps/dbsnp/*', 'steps/tools/*', '/tools/*']},

    #metadata
    author="Timothy Kuo",
    author_email="timykuo@gmail.com",
    description="Mitochondrial Genome Pipeline",
    license="MIT",
    keywords="genetics research data pipeline dependency management mitochondrial genome",
    url="https://github.com/timmykuo/mitopipeline",
    download_url="https://github.com/timmykuo/mitopipeline/releases/tag/v1.0",

    classifiers=[
        'Development Status :: 1 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
)
