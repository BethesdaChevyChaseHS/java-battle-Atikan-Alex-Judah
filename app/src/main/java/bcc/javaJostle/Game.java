package bcc.javaJostle;

import javax.swing.JPanel;
import java.awt.Graphics;
import java.io.File;
import java.lang.reflect.Constructor;
import java.net.URL;
import java.net.URLClassLoader;
import java.awt.Dimension; // For setPreferredSize
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.awt.Point; // For Point class

class Game extends JPanel {
    private ArrayList<Robot> robots;
    private ArrayList<Projectile> projectiles = new ArrayList<>(); // Placeholder
    private ArrayList<PowerUp> powerUps; // Assuming PowerUp is the class for individual power-ups
    private Map map;
    private int duration = 0;
    private Random randomGenerator; // Added for smartSpawn
    // Fields to store current display parameters
    private int currentWidth, currentHeight, currentCameraX, currentCameraY, currentTileSize;
    private double currentZoomFactor;
    private int maxDuration;

    public Game(ArrayList<String> robotFileNames, String mapName, int maxDuration) {
        this.robots = new ArrayList<>();
        this.powerUps = new ArrayList<>();
        this.maxDuration = maxDuration;
        randomGenerator = new Random(); // Initialize random generator for smartSpawn

        map = new Map(mapName);
        // use robotFileNames to create robots
       
            try {
               

                // Use the passed robotFileNames
                for (String className : robotFileNames) {
                    int encodedSpawnLocation = smartSpawn();
                    if (encodedSpawnLocation == -1) {
                        System.err.println(
                                "Could not find a valid spawn location for " + className + ". Skipping robot.");
                        continue;
                    }

                    int numCols = (this.map.getTiles() != null && this.map.getTiles().length > 0)
                            ? this.map.getTiles()[0].length
                            : 0;
                    if (numCols == 0) {
                        System.err.println("Map has no columns. Cannot place robot " + className);
                        continue;
                    }
                    int spawnRow = encodedSpawnLocation / numCols;
                    int spawnCol = encodedSpawnLocation % numCols;

                    int robotSpawnX = spawnCol * Utilities.TILE_SIZE;
                    int robotSpawnY = spawnRow * Utilities.TILE_SIZE;
                    Robot robot = Utilities.createRobot(robotSpawnX, robotSpawnY, className);
                    if (robot != null) {
                        robots.add(robot);
                        System.out.println("Added robot: " + className + " at (" + spawnCol + "," + spawnRow + ")");
                    } else {
                        System.err.println("Failed to create robot from class: " + className);
                    }

                }
                // classLoader.close(); // Consider if/when to close
            } catch (Exception e) {
                System.err.println("Failed to initialize URLClassLoader for robots directory: " + e.getMessage());
                e.printStackTrace();
            }

        // Set a preferred size for the game panel to help with layout
        setPreferredSize(new Dimension(Utilities.SCREEN_WIDTH, Utilities.SCREEN_HEIGHT));
        setFocusable(true); // Important for receiving keyboard events if needed later
    }

    // Method to update display parameters from GameManager
    public void setDisplayParameters(int width, int height, int cameraX, int cameraY, double zoomFactor) {
        this.currentWidth = width;
        this.currentHeight = height;
        this.currentCameraX = cameraX;
        this.currentCameraY = cameraY;
        this.currentZoomFactor = zoomFactor;
    }

    private int smartSpawn() {
        if (this.map == null || this.map.getTiles() == null || this.map.getTiles().length == 0
                || this.map.getTiles()[0].length == 0) {
            System.err.println("SmartSpawn: Map is not properly initialized.");
            return -1; // Invalid map
        }
        int numRows = this.map.getTiles().length;
        int numCols = this.map.getTiles()[0].length;

        List<Point> grassLocations = new ArrayList<>(); // Stores Point(col, row)
        Set<Point> visitedGrassLocations = new HashSet<>();
        int samplingAttempts = 0;
        int maxSamplingAttempts = 100; // Try up to 100 times to find 10 distinct grass spots

        while (grassLocations.size() < 10 && samplingAttempts < maxSamplingAttempts) {
            int r = randomGenerator.nextInt(numRows);
            int c = randomGenerator.nextInt(numCols);
            Point candidatePoint = new Point(c, r); // Point.x is col, Point.y is row

            if (this.map.getTiles()[r][c] == Utilities.GRASS && !visitedGrassLocations.contains(candidatePoint)) {
                grassLocations.add(candidatePoint);
                visitedGrassLocations.add(candidatePoint);
            }
            samplingAttempts++;
        }

        if (grassLocations.isEmpty()) {
            System.err.println("SmartSpawn: Could not find any grass tiles after " + maxSamplingAttempts
                    + " random attempts. Scanning map...");
            for (int r = 0; r < numRows; r++) {
                for (int c = 0; c < numCols; c++) {
                    if (this.map.getTiles()[r][c] == Utilities.GRASS) {
                        System.out.println("SmartSpawn: Found fallback grass tile at (" + c + "," + r + ")");
                        return r * numCols + c; // Return first found grass tile
                    }
                }
            }
            System.err.println("SmartSpawn: No grass tiles found on the entire map.");
            return -1; // No grass tiles anywhere
        }

        // If no robots or powerups, return the first valid grass point found
        if (this.robots.isEmpty() && this.powerUps.isEmpty()) {
            Point firstGrass = grassLocations.get(0);
            return firstGrass.y * numCols + firstGrass.x; // row * numCols + col
        }

        Point bestSpawnTile = null;
        // Initialize to a very small number, as we are looking for the maximum of these minimum distances
        double largestMinDistanceToNearestEntity = -1.0; 

        for (Point currentTile : grassLocations) {
            double tileCenterX = currentTile.x * Utilities.TILE_SIZE + Utilities.TILE_SIZE / 2.0;
            double tileCenterY = currentTile.y * Utilities.TILE_SIZE + Utilities.TILE_SIZE / 2.0;
            double currentPointMinDistanceToAnEntity = Double.MAX_VALUE;

            // Check distance to robots
            if (!this.robots.isEmpty()) {
                for (Robot robot : this.robots) {
                    if (!robot.isAlive())
                        continue;
                    double robotCenterX = robot.getX() + Utilities.ROBOT_SIZE / 2.0;
                    double robotCenterY = robot.getY() + Utilities.ROBOT_SIZE / 2.0;
                    double dist = Math.hypot(tileCenterX - robotCenterX, tileCenterY - robotCenterY);
                    currentPointMinDistanceToAnEntity = Math.min(currentPointMinDistanceToAnEntity, dist);
                }
            } else if (this.powerUps.isEmpty()) { 
                // No robots and no powerups, this case is handled earlier, but as a safeguard:
                currentPointMinDistanceToAnEntity = Double.MAX_VALUE; // Effectively makes this point desirable if it's the only one
            }


            // Check distance to powerUps
            if (!this.powerUps.isEmpty()) {
                for (PowerUp powerUp : this.powerUps) {
                    // Assuming PowerUp x,y are tile coordinates
                    double powerUpCenterX = powerUp.getX() * Utilities.TILE_SIZE + Utilities.POWER_UP_SIZE / 2.0;
                    double powerUpCenterY = powerUp.getY() * Utilities.TILE_SIZE + Utilities.POWER_UP_SIZE / 2.0;
                    double dist = Math.hypot(tileCenterX - powerUpCenterX, tileCenterY - powerUpCenterY);
                    currentPointMinDistanceToAnEntity = Math.min(currentPointMinDistanceToAnEntity, dist);
                }
            } else if (this.robots.isEmpty()) {
                 // No powerups and no robots, this case is handled earlier, but as a safeguard:
                currentPointMinDistanceToAnEntity = Double.MAX_VALUE;
            }
            
            // If there are no entities at all, currentPointMinDistanceToAnEntity will remain Double.MAX_VALUE.
            // In this scenario, any grass tile is equally good. The first one processed will be chosen.

            if (currentPointMinDistanceToAnEntity > largestMinDistanceToNearestEntity) {
                largestMinDistanceToNearestEntity = currentPointMinDistanceToAnEntity;
                bestSpawnTile = currentTile;
            }
        }

        if (bestSpawnTile != null) {
            return bestSpawnTile.y * numCols + bestSpawnTile.x; // row * numCols + col
        } else {
            // Fallback if something went wrong, though grassLocations should not be empty
            // here
            Point fallbackGrass = grassLocations.get(0);
            return fallbackGrass.y * numCols + fallbackGrass.x;
        }
    }

    @Override
    protected void paintComponent(Graphics g) {
        long startTime = System.nanoTime();
        super.paintComponent(g); // Clears the panel
        if (map != null) {
            // Pass the integer zoomFactor directly
            map.display(g, currentWidth, currentHeight, currentCameraX, currentCameraY, currentZoomFactor, robots.get(0).getDangerMap());
        }

        // Draw robots
        if (robots != null) {
            for (Robot robot : robots) {
                if (robot.getImage() != null && robot.isAlive()) {
                    int screenX = (int) (robot.getX() * currentZoomFactor - currentCameraX);
                    int screenY = (int) (robot.getY() * currentZoomFactor - currentCameraY);
                    int scaledSize = (int) (Utilities.ROBOT_SIZE * currentZoomFactor);
                    g.drawImage(robot.getImage(), screenX, screenY, scaledSize, scaledSize, null);
                }
            }
        }

        // Draw projectiles
        if (projectiles != null) {
            for (Projectile projectile : projectiles) {
                if (projectile.getProjectileImage() != null) {
                    int screenX = (int) (projectile.getX() * currentZoomFactor - currentCameraX);
                    int screenY = (int) (projectile.getY() * currentZoomFactor - currentCameraY);
                    int scaledSize = (int) (Utilities.PROJECTILE_SIZE * currentZoomFactor);
                    g.drawImage(projectile.getProjectileImage(), screenX, screenY, scaledSize, scaledSize, null);
                }
            }
        }

        // Draw powerUps
        if (powerUps != null) {
            for (PowerUp powerUp : powerUps) {
                if (powerUp.getImage() != null) {
                    double topLeftXWorld =  powerUp.getX();
                    double topLeftYWorld = powerUp.getY();

                    int screenX = (int) (topLeftXWorld * currentZoomFactor - currentCameraX);
                    int screenY = (int) (topLeftYWorld * currentZoomFactor - currentCameraY);
                    int scaledSize = (int) (Utilities.POWER_UP_SIZE * currentZoomFactor);
                    
                    g.drawImage(powerUp.getImage(), screenX, screenY, scaledSize, scaledSize, null);
                }
            }
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime); // in nanoseconds
        // System.out.println("Map drawing time: " + duration / 1_000_000.0 + " ms");
    }

    public void step() {
        final long THINK_TIME_LIMIT_MS = 5; // 5 milliseconds for robot think time

        if (robots != null) {
            for (Robot robot : robots) {
                if (!robot.isAlive()) {
                    continue; // Skip dead robots
                }

                // Execute robot.think() in a separate thread with a timeout
                Thread thinkThread = new Thread(() -> {
                    try {
                        robot.think(this.robots, projectiles, this.map, this.powerUps);
                        // If think completes without exception, it's initially considered successful
                        // The timeout check below will confirm this
                    } catch (Exception e) {
                        System.err
                                .println("Exception in Robot " + robot.getName() + " think method: " + e.getMessage());
                        e.printStackTrace();
                        robot.setSuccessfulThink(false); // Mark as unsuccessful due to exception
                    }
                });

                thinkThread.setName("RobotThinkThread-" + robot.getName());
                long startTime = System.currentTimeMillis();
                thinkThread.start();

                try {
                    thinkThread.join(THINK_TIME_LIMIT_MS); // Wait for the thread to finish with a timeout
                    long elapsedTime = System.currentTimeMillis() - startTime;

                    if (thinkThread.isAlive()) {
                        // Thread is still alive, meaning it timed out
                        thinkThread.interrupt(); // Attempt to interrupt the thread
                        // It's important that the robot's think() method checks Thread.interrupted()
                        // or handles InterruptedException to stop gracefully.
                        // We'll wait a tiny bit more for it to hopefully terminate after interrupt.
                        thinkThread.join(50); // Give it a moment to die after interrupt
                        robot.setSuccessfulThink(false);
                        System.out.println(
                                "Robot " + robot.getName() + " think method timed out after " + elapsedTime + "ms.");
                    } else if (robot.isSuccessfulThink()) {
                        // Thread finished on its own within the time limit AND no exception occurred
                        // (isSuccessfulThink would be false if an exception happened in the runnable)
                        robot.setSuccessfulThink(true);
                    }
                    // If an exception occurred in thinkThread, successfulThink is already false.
                } catch (InterruptedException e) {
                    // The Game thread itself was interrupted while waiting for the thinkThread
                    System.err.println("Game thread interrupted while waiting for robot think: " + e.getMessage());
                    robot.setSuccessfulThink(false);
                    Thread.currentThread().interrupt(); // Preserve interrupt status
                }

                // Robot.step() is final and takes Game instance
                if (robot.isAlive()) { // Robot might have been marked unsuccessful but still needs to step
                    robot.step(this);
                }
            }
        }
        //check for powerup colllisions
          if (robots != null && powerUps != null && !powerUps.isEmpty()) {
            ArrayList<PowerUp> collectedPowerUps = new ArrayList<>();
            for (Robot robot : robots) {
                if (!robot.isAlive()) {
                    continue; // Skip dead robots
                }
                // Define robot's bounding box (robot's x, y is top-left)
                double robotMinX = robot.getX();
                double robotMaxX = robot.getX() + Utilities.ROBOT_SIZE;
                double robotMinY = robot.getY();
                double robotMaxY = robot.getY() + Utilities.ROBOT_SIZE;

                for (PowerUp powerUp : powerUps) {
                    if (collectedPowerUps.contains(powerUp)) {
                        continue; // Already collected in this step by another robot
                    }

                    // Define power-up's bounding box.
                    // powerUp.getX() and powerUp.getY() are its top-left coordinates.
                    // The four corners of the power-up define its rectangular extent.
                    double p_topLeftX = powerUp.getX();
                    double p_topLeftY = powerUp.getY();
                    // Assuming POWER_UP_SIZE is for both width and height
                    double p_bottomRightX = p_topLeftX + Utilities.POWER_UP_SIZE;
                    double p_bottomRightY = p_topLeftY + Utilities.POWER_UP_SIZE;

                    // These are the min/max coordinates for the power-up's AABB
                    double powerUpMinX = p_topLeftX;
                    double powerUpMaxX = p_bottomRightX;
                    double powerUpMinY = p_topLeftY;
                    double powerUpMaxY = p_bottomRightY;

                    // Standard AABB collision check:
                    // Checks if the robot's bounding box overlaps with the power-up's bounding box.
                    // This inherently considers the full extent of the power-up defined by its corners.
                    if (robotMaxX > powerUpMinX &&    // Robot's right edge is to the right of power-up's left edge
                        robotMinX < powerUpMaxX &&    // Robot's left edge is to the left of power-up's right edge
                        robotMaxY > powerUpMinY &&    // Robot's bottom edge is below power-up's top edge
                        robotMinY < powerUpMaxY) {    // Robot's top edge is above power-up's bottom edge
                        
                        robot.applyPowerUpEffect(powerUp.getType());
                        collectedPowerUps.add(powerUp);
                        // A robot picks up only one power-up per step
                        break; 
                    }
                }
            }
            powerUps.removeAll(collectedPowerUps); // Remove all collected power-ups from the game
        }

        // Update projectiles
        if (projectiles != null) {
            for (Projectile projectile : projectiles) {
                projectile.update(this);
            }
        }

        // Remove dead projectiles
        if (projectiles != null) {
            projectiles.removeIf(projectile -> !projectile.isAlive());
        }

        // add power ups
        if(Math.random() < Utilities.POWER_UP_SPAWN_CHANCE) {
           int encodedSpawnLocation = smartSpawn(); // Use smartSpawn to find a location
           if (encodedSpawnLocation != -1) {
               int numCols = (this.map.getTiles() != null && this.map.getTiles().length > 0)
                               ? this.map.getTiles()[0].length
                               : 0;
               if (numCols > 0) {
                   int spawnRow = encodedSpawnLocation / numCols;
                   int spawnCol = encodedSpawnLocation % numCols;
                   // PowerUp constructor expects tile coordinates (col, row) or (xTile, yTile)
                   // Based on PowerUp.java, it seems x and y are tile coordinates.
                   PowerUp newPowerUp = new PowerUp((spawnCol + .5)*Utilities.TILE_SIZE -Utilities.POWER_UP_SIZE/2, (spawnRow + .5)*Utilities.TILE_SIZE - Utilities.POWER_UP_SIZE/2);
                   powerUps.add(newPowerUp);
                   System.out.println("Spawned PowerUp of type " + newPowerUp.getType() + " at tile (" + spawnCol + "," + spawnRow + ")");
               } else {
                   System.err.println("Cannot spawn power-up: Map has no columns.");
               }
           } else {
               System.err.println("Cannot spawn power-up: smartSpawn failed to find a location.");
           }
        }

        duration++;
        // System.out.println(duration + " steps taken out of " + maxDuration + " max duration.");
        
    }

    public int getDuration() {
        return duration;
    }

    public int getMaxDuration() {
        return maxDuration;
    }

    public ArrayList<Robot> getRobots() {
        return robots;
    }

    public void addProjectile(Projectile projectile) {
        projectiles.add(projectile);
    }

    public Map getMap() {
        return map;
    }

    public ArrayList<Projectile> getProjectiles() {
        return projectiles;
    }

    public ArrayList<PowerUp> getPowerUps() {
        return powerUps;
    }

    public boolean isGameOver() {
        // Check if the game duration has reached the maximum allowed duration
        if (duration >= maxDuration) {
            return true; // Game is over due to time limit
        }

        // Check if exactly one robot is alive
        int aliveCount = 0;
        for (Robot robot : robots) {
            if (robot.isAlive()) {
                aliveCount++;
            }
        }
        return aliveCount == 1;
    }

    public Robot getWinner() {
        // If the game is over, return the only remaining alive robot
        if (!isGameOver()) {
            return null; // Game is not over yet
        }
        if (duration >= maxDuration) {
            int highestHealthPercentage = -1;
            Robot winner = null;
            for (Robot robot : robots) {
                if (robot.isAlive()) {
                    int healthPercentage = (int) ((robot.getHealth() / (double) robot.getMaxHealth()) * 100);
                    if (healthPercentage > highestHealthPercentage) {
                        highestHealthPercentage = healthPercentage;
                        winner = robot;
                    }
                }
            }
            return winner; // Return the robot with the highest health percentage
        }
        for (Robot robot : robots) {
            if (robot.isAlive()) {
                return robot; // Return the winning robot
            }
        }
        return null; // No winner found, should not happen if game is over

    }
}