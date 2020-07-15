from setuptools import setup, find_packages
setup(
    name="dao_mop",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "dao_mop_sns = dao_mop.sns:main",
            "dao_mop_train_cnn = dao_mop.train_model:main",
            "dao_mop_build_plant_db = dao_mop.build_plant_list_db:main",
        ],
    }
)
