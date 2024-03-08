

Referencias:
1. https://carpentries-incubator.github.io/workflows-nextflow/01-getting-started-with-nextflow.html


- What is a workflow and what are workflow management systems?

Analysing data involves a sequence of tasks, including gathering, cleaning, and processing data. *Workflow Management Systems (WfMS)* such as Snakemake, Galaxy, and Nextflow have been developed specifically to manage computational data-analysis workflows. 

- Why should I use a workflow management system?

Provide useful features to improve

Run time management; Software management; Portability & Interoperability; Reproducibility; Re-entrancy

- What are the main features of Nextflow?

*Fast prototyping:* reuse existing scripts and tools;

*Reproducibility:* Supports several container technologies, such as Docker and  Singularity, as well as the package manager Conda.

*Portability and Interoperability*: separates functional logic (the steps of the workflow) from the execution settings (how the workflow is executed). THis allows the pipeline to be run on multiple platforms (without changing the steps of the workflow)

*Simple parallelism*: greatly simplifies the splitting of tasks that can be run at the same time (parallelisation)

*Continuous checkpoints & re-entrancy*: all steps tracked

- What are the main components of a Nextflow script?

*Processes:* describe a task to be run

*Channels:* used to manipulate the flow of data from one process to the next

*Workflows:* define the interaction between process and ultimately the pipeline execution flow itself.

![Alt text](<Screenshot from 2023-12-28 22-08-55.png>) [1]


### TODO

- createa yaml/json file to input parameters (best practices for sveral values)
-