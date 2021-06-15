from setuptools import setup, find_packages

# install the daomop! -- test

version = {}
with open("daomop/version.py") as fp:
    exec(fp.read(), version)
print(version['__version__'])

install_requires=[
'numpy~=1.20.1',
'astropy~=4.2',
'pyraf~=2.1.15',
'vos~=3.3',
'matplotlib~=3.3.4',
'scipy~=1.6.1',
'ephem~=3.7.7.1',
'pyds9~=1.8.1',
'daomop~=0.3.3',
'setuptools~=56.0.0',
'ccdproc~=2.1.1']

setup(
    name="daomop",
    version=version['__version__'],
    packages=find_packages(),
    scripts=['scripts/daomop-stack-cmd',
             'scripts/daomop-run-sns.sh',
             'scripts/daomop-target-sns.sh',
             'scripts/daomop-get-diffs.sh',
             'scripts/daomop-maskdiffs.sh',
             'scripts/daomop-filesync.sh'],
    package_data={'daomop': ['*.json']},
    entry_points={
        "console_scripts": [
            "daomop-sns = daomop.sns:main",
            "daomop-predict-obs = daomop.predict_obs:run",
            "daomop-measure = daomop.measure_kbo:run",
            "daomop-measure-orbit = daomop.measure_from_orbit:run",
            "daomop-train-cnn = daomop.train_model:main",
            "daomop-build-plant-db = daomop.build_plant_list_db:main",
            "daomop-intelligentMasker = daomop.intelligentMasker:main",
        ],
    }
)
