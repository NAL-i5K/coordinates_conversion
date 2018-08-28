from setuptools import setup, find_packages

setup(
    name='coordinates_conversion',
    author='NAL i5k workspace',
    author_email='i5k@ars.usda.gov',
    packages=find_packages(),
    install_requires=['pysam'],
    entry_points={
        'console_scripts': [
            'fasta_diff=coordinates_conversion.fasta_diff:main',
            'update_bam=coordinates_conversion.update_bam:main',
            'update_bed=coordinates_conversion.update_bed:main',
            'update_bedgraph=coordinates_conversion.update_bedgraph:main',
            'update_gff=coordinates_conversion.update_gff:main',
            'update_vcf=coordinates_conversion.update_vcf:main',
        ]
    }
)
