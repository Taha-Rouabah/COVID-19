
import setuptools


def get_requirements():
    reqs = []
    for line in open('requirements.txt', 'r').readlines():
        reqs.append(line)
    return reqs


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="seiqrdp_model",
    version="8.0-alpha",
    author="T.Rouabah, N.Belaloui, A.Tounsi",
    author_email="m.t.rouabah@gmail.com",
    description="Computational tools to fit data, simulate and calibrate\
    parameters of the compartmental epidemiological SEIQRDP model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Taha-Rouabah/COVID-19",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=get_requirements()
)
