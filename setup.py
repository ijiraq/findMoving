from setuptools import setup, find_packages

version = {}
with open("daomop/version.py") as fp:
    exec(fp.read(), version)
print(version['__version__'])

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
