from setuptools import setup, find_packages
setup(
    name="daomop",
    version="0.1",
    packages=find_packages(),
    scripts=['scripts/daomop-stack-cmd'],
    entry_points={
        "console_scripts": [
            "daomop-sns = daomop.sns:main",
            "daomop-train-cnn = daomop.train_model:main",
            "daomop-build-plant-db = daomop.build_plant_list_db:main",
        ],
    }
)
