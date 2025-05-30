[![Open in Codespaces](https://classroom.github.com/assets/launch-codespace-2972f46106e565e64193e422d61a12cf1da4916b45550586e14ef0a7c637dd04.svg)](https://classroom.github.com/open-in-codespaces?assignment_repo_id=19673592)
# Setup Instructions
In order to run the GUI:
Go to the ports tab. Click "Add port" and type the number 6080.
Follow the link to your forwarded address. This is where your GUI should be available once you run your code. Note that the GUI won't work very well until you have finished most of the code - tests will be used to make sure the intermediate steps works well. 

To access the code, go to app > src > main > java


# Programming Instructions

General notes:
1. Don't change any method headers. These will be correct. 
2. For testing run the given command in the terminal. It will be of the pattern gradle ______
3. Say hello to your (maybe) first real code editor. Autocomplete, error highlighting, and other things are about to make your life so much easier. 
4. Apply for the student github developer pack!!


# Submission
1. Testing everything. Run the command - 
``` 
gradle test
```

2. Submit on github classroom by running the following commands. This also saves your work permanently (unless you actively want to delete it). 

```
git add . 
git commit -m "submitting"
git push
```

# DELETE THIS: 
1. Clone existing repo with desktoplite
2. Run `gradle init --type java-application --package bcc.projectname`
3. Select version 11 for Java (CSA version), select groovy and JUnit 4
4. Add extension for java
5. Add the following to build.gradle for verbose testing: 
```
test {
    testLogging {
        events "passed", "failed", "skipped"
        exceptionFormat "full"
        showStandardStreams = true
    }
}

tasks.withType(Test) {
    testLogging {
        events "passed", "failed", "skipped"
        exceptionFormat "full"
        showStandardStreams = true
    }
}
```
6. For additional test files (recommended for more granular grading):
```
task runTestName(type: Test) {
    description = "Runs tests in the testName class."
    group = "verification"
    filter {
        includeTestsMatching "bcc.packagename.testName"
    }
}
```