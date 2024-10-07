import pandas as pd

# Create the new schedule for September 2024 using the task rotation among new roommates
dates = ["Oct 1-6", "Oct 7-13", "Oct 14-20", "Oct 21-27", "Oct 28-31"]
tasks = [
    ["Kitchen cleaning", "Trash curbside pickup"],
    ["Hall, living room Vacuuming", "Trash curbside pickug"],
    ["Kitchen cleaning", "Trash curbside pickup"],
    ["Hall, living room Vacuuming", "Trash curbside pickug"],
    ["Kitchen cleaning", "Trash curbside pickup"s]
]

# Roommates cycling through the tasks
roommates = ["Senthil", "June", "Nicole", "Chris Ruping", "Pranav"]

# Creating DataFrame
new_schedule = pd.DataFrame({
    "Date": dates,
    "Task 1": [task[0] for task in tasks],
    "Roommate 1": [roommates[i % len(roommates)] for i in range(len(dates))],
    "Task 2": [task[1] for task in tasks],
    "Roommate 2": [roommates[(i + 1) % len(roommates)] for i in range(len(dates))]
})

new_schedule_path = '/Users/sq566/Desktop/New_House_Cleaning_Schedule_October_2024.xlsx'
new_schedule.to_excel(new_schedule_path, index=False)

new_schedule, new_schedule_path
