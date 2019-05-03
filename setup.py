from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="mitopipeline",
    version="1.0.post1",
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
    long_description=long_description,
    long_description_content_type="text/x-rst",
    keywords="genetics research data pipeline dependency management mitochondrial genome",
    url="https://github.com/timmykuo/mitopipeline",
    download_url="https://github.com/timmykuo/mitopipeline/releases/tag/v1.0",

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Sociology :: Genealogy'
    ],
)
