# assume we are a subdir of the same base dir as 'Teaching' repository
import os
import shutil
import glob

courses = {
	"MASx50" : {
		"notes.pdf" : "notes",
	},
	"MASx52" : {
		"formula_sheet.pdf" : os.path.join("notes","formula_sheet"),
		"notes_0.pdf" : os.path.join("notes","notes_0"),
		"notes_1.pdf" : "notes",
		"notes_2.pdf" : "notes",
	}
}

for course, material in courses.items():
		for name, folder in material.items():
			in_path = os.path.join("..", "Teaching", course, folder, name)
			out_path = os.path.join(".", course, name)
			print(in_path, out_path)
			shutil.copy(in_path, out_path)
			
		in_html_folder = os.path.join("..", "Teaching", course, "notes", "html_"+course)
		out_html_folder = os.path.join(".", course, "html")
		if os.path.exists(out_html_folder):
			shutil.rmtree(out_html_folder)
		print(in_html_folder, out_html_folder)
		shutil.copytree(in_html_folder, out_html_folder)