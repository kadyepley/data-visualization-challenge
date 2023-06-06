# data-visualization-challenge# Pymaceuticals Inc.
---

### Analysis

- Add your analysis here.
 
# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
combined_data =  pd.merge(mouse_metadata, study_results,how='right') # choosing 'right' to use the keys from mouse_metadata

# Display the data table for preview
combined_data = combined_data[['Mouse ID', 'Timepoint', 'Tumor Volume (mm3)', 'Metastatic Sites', 'Drug Regimen', 'Sex', 'Age_months', 'Weight (g)']] # We use the double bracket so the output is a dataframe

combined_data.head ()
# Checking the number of mice.
mice = combined_data["Mouse ID"].value_counts() # Assigning a variable to count under Mouse ID
mice_count = len(mice) # Actually counting all the Mouse IDs
mice_count # Showing the outcome of 249 mouse IDs
# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
dup_id = combined_data.loc[combined_data.duplicated(subset=['Mouse ID','Timepoint',]),'Mouse ID'].unique() # Access the dataframe and search for duplicates of the Mouse IDs
dup_id # Showing the outcome of this search
# Optional: Get all the data for the duplicate mouse ID. 
dup_id_df = combined_data.loc[combined_data["Mouse ID"] == "g989", :] # Setting dup_id_df to equal all the lines (and their data) for each instance of Mouse ID = g989
dup_id_df # Showing this dataframe
# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = combined_data[combined_data['Mouse ID'].isin(dup_id) == False] # Creating df with all instances where the Mouse ID in the combined data frame is not equal to the duplicated ID 
clean_df.head () # Showing this dataframe
# Checking the number of mice in the clean DataFrame.
mice_clean = clean_df["Mouse ID"].value_counts() # Assigning a variable to count the Mouse IDs in the clean_df
mice_clean_count = len(mice_clean) # Actually counting all the Mouse IDs
mice_clean_count # Showing the outcome of 248 mouse IDs. This is just one less than the previous count because the dup_id was removed

## Summary Statistics
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.

drug_reg_mean = clean_df.groupby('Drug Regimen').mean()["Tumor Volume (mm3)"]
drug_reg_med = clean_df.groupby('Drug Regimen').median()["Tumor Volume (mm3)"]
drug_reg_var = clean_df.groupby('Drug Regimen').var()["Tumor Volume (mm3)"]
drug_reg_stdev = clean_df.groupby('Drug Regimen').std()["Tumor Volume (mm3)"]
drug_reg_sem = clean_df.groupby('Drug Regimen').sem()["Tumor Volume (mm3)"]

# Create each grouping of mean, median, variance, std dev, and std error for each drug name to be pulled in to a dataframe

summary_stat_table = pd.DataFrame({"Mean Tumor Volume": drug_reg_mean,"Median Tumor Volume": drug_reg_med,"Tumor Volume Variance": drug_reg_var,"Tumor Volume Std. Dev.": drug_reg_stdev,"Tumor Volume Std. Err.": drug_reg_sem})

summary_stat_table
# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line

summary_agg_table = clean_df.groupby(['Drug Regimen'])[['Tumor Volume (mm3)']].agg(['mean','median', 'var', 'std', 'sem']) 
summary_agg_table
## Bar and Pie Charts
# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.

mice_count = clean_df["Drug Regimen"].value_counts() # Count the mice using each drug 
plot_pandas = mice_count.plot.bar(color = 'tab:blue') # Choosing the bar graph with blue bars
plt.xlabel("Drug Regimen") # Name the x-axis
plt.ylabel("Number of Observed Mouse Timepoints") # Name the y-axis

# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.

x_axis = mice_count.index.values
y_axis = mice_count.values

plt.bar(x_axis, y_axis, color = 'tab:blue', alpha = 0.8, align = 'center')

plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timepoints")
plt.xticks(rotation = "vertical")
# Generate a pie plot showing the distribution of female versus male mice using Pandas
data_by_sex = clean_df["Sex"].value_counts() # Counting the number of each sex in the cleaned data
data_by_sex.plot.pie(autopct = "%1.1f%%") # Plotting as a pie chart with one decimal point in the percentages
# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ['Female', 'Male'] # Set the labels
sizes = [49.0, 51.0] # Manually setting the sizes of each side of the pie chart
plot = data_by_sex.plot.pie(y = 'Total Count', autopct = "%1.1f%%")
plt.ylabel('Sex')
## Quartiles, Outliers and Boxplots
# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
each_mouse = combined_data.groupby(["Mouse ID"]).max()
each_mouse_cleared = each_mouse.reset_index()

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_df = each_mouse_cleared[['Mouse ID', 'Timepoint']].\
    merge(combined_data, on = ['Mouse ID', 'Timepoint'], how = "left")



# Put treatments into a list for for loop (and later for plot labels)
#treatments = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
#drugs = combined_data[combined_data["Drug Regimen"].isin(treatments)]

# Calculate the IQR and quantitatively determine if there are any potential outliers. 

# Not going to lie, trying to get this to match the outcomes and have the prompts met confused me a bit. 
def get_outliers(regimen):
    regimen_df = merged_df.loc[merged_df["Drug Regimen"] == regimen]['Tumor Volume (mm3)']
    
    quartiles = regimen_df.quantile([.25,.5,.75])
    quart_first = quartiles[0.25]
    quart_last = quartiles[0.75]
    quart_range = quart_last - quart_first
    lower_bound = quart_first - (1.5 * quart_range)
    upper_bound = quart_last + (1.5 * quart_range)

    outliers = regimen_df.loc[(regimen_df < lower_bound) | (regimen_df > upper_bound)]

    print(f"{regimen}'s potential outliers: {outliers}")
    return regimen_df

capo = get_outliers("Capomulin")
rami = get_outliers("Ramicane")
infu = get_outliers("Infubinol")
ceft = get_outliers("Ceftamin")
    
    # Locate the rows which contain mice on each drug and get the tumor volumes

    # add subset 
    
    # Determine outliers using upper and lower bounds

# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
red_dot = dict(markerfacecolor = "red")
plt.boxplot([capo, rami, infu, ceft], labels = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin'], flierprops = red_dot)
plt.ylabel('Final Tumor Volume (mm3)')
plt.show()
## Line and Scatter Plots
# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_df = clean_df.loc[clean_df["Drug Regimen"] == "Capomulin",:] 
line_df = capomulin_df.loc[capomulin_df["Mouse ID"] == "l509",:]
line_df.head()
x_axis = line_df["Timepoint"]
tumor_size = line_df["Tumor Volume (mm3)"]
fig1, ax1 = plt.subplots()

# Format and labels

plt.title('Capomulin treatment of moust l509')
plt.plot(x_axis, tumor_size, linewidth = 2, color = "steelblue", label = "Fahreneit")
plt.xlabel('Timepoint (days)')
plt.ylabel('Tumor Volume (mm3)')
plt.show()

# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
fig1, ax1 = plt.subplots()
avg_capo_vol = capomulin_df.groupby(['Mouse ID']).mean()


#Format
plt.scatter(avg_capo_vol['Weight (g)'], avg_capo_vol['Tumor Volume (mm3)'], color = "steelblue")
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
marker_size = 15
plt.show()
## Correlation and Regression
# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
corre_coeff = st.pearsonr(avg_capo_vol['Weight (g)'], avg_capo_vol['Tumor Volume (mm3)'])
print(f"The correlation between mouse weight and the average tumor volume is {round(corre_coeff[0], 2)}")

# Calculating the line of best fit
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(avg_capo_vol["Weight (g)"], avg_capo_vol["Tumor Volume (mm3)"])
values_for_formula = avg_capo_vol["Weight (g)"] * slope + intercept

plt.scatter(avg_capo_vol["Weight (g)"], avg_capo_vol["Tumor Volume (mm3)"], color = 'steelblue')
plt.plot(avg_capo_vol["Weight (g)"], values_for_formula, color = 'red')
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.show()
