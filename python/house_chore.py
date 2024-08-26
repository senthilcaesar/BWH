import pandas as pd

# Create the new schedule for September 2024 using the task rotation among new roommates
dates = ["Sept 1-7", "Sept 8-14", "Sept 15-21", "Sept 22-28", "Sept 29-30"]
tasks = [
    ["Kitchen cleaning", "Trash curbside pickup"],
    ["Hall, living room Vacuuming", "stairs, dust cleaning"],
    ["Kitchen cleaning", "Trash curbside pickup"],
    ["Hall, living room Vacuuming", "stairs, dust cleaning"],
    ["Kitchen cleaning", "Trash curbside pickup"]
]

# Roommates cycling through the tasks
roommates = ["Senthil", "June", "Nicole", "Chris", "Pranav"]

# Creating DataFrame
new_schedule = pd.DataFrame({
    "Date": dates,
    "Task 1": [task[0] for task in tasks],
    "Roommate 1": [roommates[i % len(roommates)] for i in range(len(dates))],
    "Task 2": [task[1] for task in tasks],
    "Roommate 2": [roommates[(i + 1) % len(roommates)] for i in range(len(dates))]
})

new_schedule_path = '/Users/sq566/Desktop/New_House_Cleaning_Schedule_September_2024.xlsx'
new_schedule.to_excel(new_schedule_path, index=False)

new_schedule, new_schedule_path
