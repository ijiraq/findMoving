from setuptools import setup, find_packages
setup(
    name="dao_mop",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "shift_and_stack = dao_mop.sns:main",
            "train_mop_model = dao_mop.train_model:main",
            "build_plant_db = dao_mop.build_plant_list_db:main",
        ],
    }
)
