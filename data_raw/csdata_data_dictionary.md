Dataset is a sample of tow-away crashes for vehicles less than 10 years old in the
united states. Sample is constructed by selecting PSU, then selecting police jurisdiction
within PSU, then serious injuries are oversampled. So the data must be re-weighted
relative to ratwgt to be representative of the population of accidents in the US.

# Notes on data usability/quirks

* There seems to be some errors on height and weight, since I calculate BMI over
  85, the max any table shows.

# Variable identification
## Identifier variables
Taken together, these variables uniquely identify each row

* gyear - Year of the accident
* psu - Primary sampling unit - 27 sites around the country
* caseno - Identifier for the accident
* vehno - Identifier for the vehicle in an accident
* occno - Identifier for an occupant of a vehicle in an accident

## Other "administrative" variables
* ratwgt - Case weight. Some accidents are over-sampled. Some under-sampled.
           Number of cases in a given year that a single sampled case represents. 
           Repeating 

## Crash characteristic variables
* mais - Maximum Abbreviated Injury Scale Level, severity of injury. mais >= 3 is severe
* bmode - Crash mode. Frontal, Rear, far side (FS), near side (NS), roll, unknown
          FS and NS depend on occupant seating. For example, a crash to the passenger
          side is considered FS for occupants on the drivers side and NS for occupants 
          on the passenger side.
* dvtotal - measure of severity of a traffic collision, defined as the change
            in velocity between pre-collision and post-collision trajectories of 
            a vehicle.
* delta_v - scaled measure of dvtotal. 
* maxc1 - maximum exterior crush for the first impact area in inches
* maxc2 - maximum exterior crush for the second impact area in inches
* intrus_p - maximum depth of intrusion near a given occupant defined by occno
* multi - flag for multiple vehicle collison. 0 - No, 1 - Yes
* rollover - Number of 1/4 rollover turns that a vehicle makes
* cdc_typ - Collision Deformation Classification (categorical)
* pdo1 - Principal Direction of Force
* shl1 - Specific horizontal location (categorical)
* dvl - Width of damage
* dvd - Offset of damage from the center lines of the vehicle

## Occupant characteristic variables
* belt_use - Indicator for passenger seatbelt use
* age - Age of occupant
* height - height of occupant in centimeters
* weight - weight of occupant in Kg
* bmi - Body Mass Index: weight / height ^ 2 
* sex - Categorical variable 
* seatrack - Seat track position
* seatpos - seating position

## Vehicle characteristic variables
* curbwgt_kg - weight of the vehicle in Kg.
* btype - Body type of the vehicle: 1 - Car, 2 - SUV, 3 - Pickup Truck, 4 - Van
* makmod - Numeric identifier for vehicle make and model
* modelyr - Model year of vehicle

# Additional data notes to go in paper

## Data Source
(http://www.nhtsa.giov/Data/National+Automotive+Sampling+System+(NASS)/NASS+Crashworthiness+Data+System)

For this study, data are compiled from the National Automotive Sampling System / Crashworthiness Data System (NASS/CDS). NASS/CDS data are collected by the National Highway Traffic Safety Administration. The data represent a probability sample of crashes in the United States involving passenger cars, light trucks, vans, and utility vehicles. Approximately 5,000 accidents are sampled per year, stratified by primary sampling unit (PSU). There are 27 PSUs across the country. Each PSU has an assigned crash investigation team. These crash investigators obtain data by investigating the crash site, studying the vehicles involved, interviewing crash victims, by reviewing medical records. Personal information about individuals involved in the crashes are not included in any public NASS file.

The data selected for this study include NASS/CDS years 2000 - 2014. Several observations are removed due to data misingess. [Need more information here.]

## Sample Weighting

NASS/CDS data is a probability sample of crashes in the United States, but it is not representative. More severe crashes are prioritized for inclusion in the database. Therefore each accident is given a sample weight, proportional to the number of crashes in the population it represents. 

Instead of analyzing the raw NASS/CDS data, this study constructs two analysis datasets. The first is a resampling of the data to be representative of the population of crashes in the United States. A bootstrapped sample is taken where each observation's probability of being sampled is proportional to its sampling weight. Some cases have abnormally large sampling weights, however. These tend to be cases that are very rare in the data set, but very common in a particular geographic region. To avoid these singular cases being oversampled, with a few cases dominating the data set, any weight larger than [3,400] is recoded to be [3,400] before bootstrapping. The second analysis dataset is used to study severe injuries. A bootstrapped sample is where each observation's probability of being sampled is proportional to MAIS, excluding observations where MAIS is unknown (coded as "7" or missing).

## Maximum Abbreviated Injury Scale (MAIS)

(Source: Baker SP, O'Neill B, Haddon W, Long WB, "The Injury Severity Score: a method for describing patients with multiple injuries and evaluating emergency care,", Journal of Trauma, Vol. 14, No. 5, 1997, pp. 187-196.)

The Abbreviated Injury Scale (AIS) is a coding system created by the Association for the Advancement of Automotive Medicine. It is used to classify and describe  the severity of injuries, representing the threat to life associated with an injury, not a comprehensive assesement of injury severity. The AIS scale is listed in the table below. 

| AIS-Code | Injury Level | AIS % Prob. of Death |
|----------|--------------|----------------------|
| 1        | Minor        | 0                    |
| 2        | Moderate     | 1 - 2                |
| 3        | Serious      | 8 - 10               |
| 4        | Severe       | 5 - 50               |
| 5        | Critical     | 5 - 50               |
| 6        | Maximum      | 100                  |
| 9        | Not specified| N/A                  |


In an automobile accident, occupants may have several injuries. Therefore, the maximum AIS (MAIS) is recorded as the injury severity for an occupant in our data. For example, if an occupant has a superficial laceration of the thigh (AIS level 1) and a fractured sternum (AIS level 2), the MAIS is 2. Additionally, the MAIS  variable in the NASS/CDS data is coded as "7" where AIS would be "9". A common proxy for severe injuries in the literature is "MAIS 3+", a binary indicator for whether an occupant's MAIS injury is in the range 3 - 6. 

