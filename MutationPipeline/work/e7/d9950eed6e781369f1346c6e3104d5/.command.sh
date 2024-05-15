#!/bin/bash -ue
# Run the getvariants.R script using the provided renviron file, start-to-end date and city name

getvariants.R --renviron=.Renviron --start_collection=2020-01-01 --end_collection=2023-01-01 --location='South America / Brazil / Parana / Foz'
