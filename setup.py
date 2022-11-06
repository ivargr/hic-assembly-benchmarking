from setuptools import setup

setup(name='hic_phaser',
      version='0.0.1',
      description='hic_phaser',
      url='http://github.com/ivargr/hic_phaser',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=["hic_phaser"],
      zip_safe=False,
      install_requires=['numpy', 'shared_memory_wrapper>=0.0.8', 'bionumpy==0.2.2', 'typer'],
      classifiers=[
          'Programming Language :: Python :: 3'
      ],
      entry_points={
          'console_scripts': ['hic_phaser=hic_phaser.command_line_interface:main']
      }
      )
