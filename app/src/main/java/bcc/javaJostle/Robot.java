package bcc.javaJostle;

import java.util.ArrayList;
import java.awt.image.BufferedImage;

public abstract class Robot {
    // attribute points
    private int healthPoints;
    private int speedPoints;
    private int attackSpeedPoints;
    private int projectileStrengthPoints;

    // attribute calculated values
    private int health;
    private int maxHealth;
    private int speed; // This is an int
    private int attackMaxCooldown;
    protected int attackCurCooldown;
    private int projectileSpeed;
    private int projectileDamage;
    private String name;
    private BufferedImage image;
    private BufferedImage projectileImage;
    private boolean successfulThink = true;

    private int x;
    private int y;

    protected int xMovement;
    protected int yMovement;
    protected int xTarget;
    protected int yTarget;
    protected boolean shoot;

    // Power-up effect fields
    private static final int BOOST_DURATION_TICKS = 150;
    private int speedBoostDuration = 0;
    private int attackBoostDuration = 0;

    private int originalSpeed;
    private int originalProjectileSpeed;
    private int originalProjectileDamage;


    public Robot(int x, int y, int healthPoints, int speedPoints, int attackSpeedPoints, int projectileStrengthPoints,
            String robotName, String imageName, String projectileImageName) {
        // all values need to be from 1-5, summing to 10 in total
        int sum = healthPoints + speedPoints + attackSpeedPoints + projectileStrengthPoints;

        if (healthPoints < 1) {
            throw new IllegalArgumentException("Health points must be at least 1");
        } else if (speedPoints < 1) {
            throw new IllegalArgumentException("Speed points must be at least 1");
        } else if (attackSpeedPoints < 1) {
            throw new IllegalArgumentException("Attack speed points must be at least 1");
        } else if (projectileStrengthPoints < 1) {
            throw new IllegalArgumentException("Projectile strength points must be at least 1");
        } else if (sum != 10) {
            throw new IllegalArgumentException("The sum of all points must equal 11");
        } else if (healthPoints > 5) {
            throw new IllegalArgumentException("Health points must not exceed 5");
        } else if (speedPoints > 5) {
            throw new IllegalArgumentException("Speed points must not exceed 5");
        } else if (attackSpeedPoints > 5) {
            throw new IllegalArgumentException("Attack speed points must not exceed 5");
        } else if (projectileStrengthPoints > 5) {
            throw new IllegalArgumentException("Projectile strength points must not exceed 5");
        }
        try {
            this.image = Utilities.loadImage(imageName);
            if (this.image == null)
                this.image = Utilities.ROBOT_ERROR;
        } catch (Exception e) {
            System.err.println("Failed to load image: " + imageName + " - " + e.getMessage());
            e.printStackTrace();
            this.image = Utilities.ROBOT_ERROR;
        }
        try {
            this.projectileImage = Utilities.loadImage(projectileImageName);
            if (this.projectileImage == null)
                this.projectileImage = Utilities.DEFAULT_PROJECTILE_IMAGE;
        } catch (Exception e) {
            System.err.println("Failed to load image: " + projectileImageName + " - " + e.getMessage());
            e.printStackTrace();
            this.projectileImage = Utilities.DEFAULT_PROJECTILE_IMAGE;
        }

        this.x = x;
        this.y = y;
        this.name = robotName;

        this.healthPoints = healthPoints;
        this.speedPoints = speedPoints;
        this.attackSpeedPoints = attackSpeedPoints;
        this.projectileStrengthPoints = projectileStrengthPoints;

        // give stats based on points
        this.health = 30 + healthPoints * 20; // 50 - 130
        this.maxHealth = this.health;
        this.speed = 2 + speedPoints; // 3 - 7 (int)
        this.attackMaxCooldown = 22 - attackSpeedPoints * 2; // 20 - 12
        this.attackCurCooldown = attackMaxCooldown;
        this.projectileSpeed = 5 + projectileStrengthPoints; // 6 - 10
        this.projectileDamage = 10 + projectileStrengthPoints * 3; // 13 - 25

        this.originalSpeed = this.speed;
        this.originalProjectileSpeed = this.projectileSpeed;
        this.originalProjectileDamage = this.projectileDamage;
    }

    public boolean canAttack() {
        return attackCurCooldown <= 0;
    }

    protected void shootAtLocation(int x, int y) {
        xTarget = x;
        yTarget = y;
        shoot = true;
    }

    public boolean[][] getDangerMap() {
        return new boolean[1][1];
    }

    public abstract void think(ArrayList<Robot> robots, ArrayList<Projectile> projectiles, Map map,
            ArrayList<PowerUp> powerups);

    private boolean isPointOkay(int pX, int pY, Map gameMap, ArrayList<Robot> allRobots) {
        if (gameMap == null || gameMap.getTiles() == null)
            return false;
        int[][] mapTiles = gameMap.getTiles();
        int mapRows = mapTiles.length;
        if (mapRows == 0)
            return false;
        int mapCols = mapTiles[0].length;
        if (mapCols == 0)
            return false;

        // 1. Map Boundary Check for the point
        if (pX < 0 || pY < 0 || pX >= mapCols * Utilities.TILE_SIZE || pY >= mapRows * Utilities.TILE_SIZE) {
            return false;
        }

        // 2. Wall Tile Check for the point
        int tileCol = (int) (pX / Utilities.TILE_SIZE);
        int tileRow = (int) (pY / Utilities.TILE_SIZE);
        if (tileRow < 0 || tileRow >= mapRows || tileCol < 0 || tileCol >= mapCols) {
            return false;
        }
        if (mapTiles[tileRow][tileCol] == Utilities.WALL) {
            return false;
        }

        // 3. Other Robot Check for the point
        if (allRobots != null) {
            for (Robot otherRobot : allRobots) {
                if (otherRobot == this || !otherRobot.isAlive()) {
                    continue;
                }
                if (pX >= otherRobot.getX() && pX < (otherRobot.getX() + Utilities.ROBOT_SIZE) &&
                        pY >= otherRobot.getY() && pY < (otherRobot.getY() + Utilities.ROBOT_SIZE)) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean canMoveTo(int targetX, int targetY, Game game) {
        // Define the 4 corners at this potential new position
        // Top-left, Top-right, Bottom-left, Bottom-right
        int c1x = targetX;
        int c1y = targetY;
        int c2x = targetX + Utilities.ROBOT_SIZE - 1; // Use -1 for inclusive edge if ROBOT_SIZE is a dimension
        int c2y = targetY;
        int c3x = targetX;
        int c3y = targetY + Utilities.ROBOT_SIZE - 1;
        int c4x = targetX + Utilities.ROBOT_SIZE - 1;
        int c4y = targetY + Utilities.ROBOT_SIZE - 1;

        return isPointOkay(c1x, c1y, game.getMap(), game.getRobots()) &&
                isPointOkay(c2x, c2y, game.getMap(), game.getRobots()) &&
                isPointOkay(c3x, c3y, game.getMap(), game.getRobots()) &&
                isPointOkay(c4x, c4y, game.getMap(), game.getRobots());
    }

    private boolean isTileMud(int currentX, int currentY, Map gameMap) {
        if (gameMap == null || gameMap.getTiles() == null)
            return false;
        int[][] mapTiles = gameMap.getTiles();
        int mapRows = mapTiles.length;
        if (mapRows == 0)
            return false;
        int mapCols = mapTiles[0].length;
        if (mapCols == 0)
            return false;

        int[] cornersX = { currentX, currentX + Utilities.ROBOT_SIZE - 1, currentX,
                currentX + Utilities.ROBOT_SIZE - 1 };
        int[] cornersY = { currentY, currentY, currentY + Utilities.ROBOT_SIZE - 1,
                currentY + Utilities.ROBOT_SIZE - 1 };

        for (int i = 0; i < 4; i++) {
            int tileCol = (int) (cornersX[i] / Utilities.TILE_SIZE);
            int tileRow = (int) (cornersY[i] / Utilities.TILE_SIZE);

            if (tileRow >= 0 && tileRow < mapRows && tileCol >= 0 && tileCol < mapCols) {
                if (mapTiles[tileRow][tileCol] == Utilities.MUD) {
                    return true;
                }
            }
        }
        return false;
    }
    public void applyPowerUpEffect(String type) {
        System.out.println(this.name + " picked up " + type + " power-up!");
        switch (type) {
            case "health":
                this.health += this.maxHealth / 2;
                System.out.println(this.name + " new health: " + this.health);
                break;
            case "speed":
                if (speedBoostDuration == 0) { // Only apply if not already boosted
                    this.originalSpeed = this.speed; // Store current speed if it was somehow changed by other means
                }
                this.speed = this.originalSpeed * 2; // Apply boost based on original
                this.speedBoostDuration = BOOST_DURATION_TICKS;
                System.out.println(this.name + " new speed: " + this.speed + " for " + BOOST_DURATION_TICKS + " ticks.");
                break;
            case "attack":
                if (attackBoostDuration == 0) { // Only apply if not already boosted
                    this.originalProjectileSpeed = this.projectileSpeed;
                    this.originalProjectileDamage = this.projectileDamage;
                }
                this.projectileSpeed = this.originalProjectileSpeed * 2;
                this.projectileDamage = this.originalProjectileDamage * 2;
                this.attackBoostDuration = BOOST_DURATION_TICKS;
                System.out.println(this.name + " new projectile speed: " + this.projectileSpeed + ", new damage: " + this.projectileDamage + " for " + BOOST_DURATION_TICKS + " ticks.");
                break;
            default:
                System.err.println("Unknown power-up type: " + type);
                break;
        }
    }

    private void updatePowerUpEffects() {
        if (speedBoostDuration > 0) {
            speedBoostDuration--;
            if (speedBoostDuration == 0) {
                this.speed = this.originalSpeed; // Revert to original speed
                System.out.println(this.name + " speed boost wore off. Speed reverted to " + this.speed);
            }
        }
        if (attackBoostDuration > 0) {
            attackBoostDuration--;
            if (attackBoostDuration == 0) {
                this.projectileSpeed = this.originalProjectileSpeed; // Revert
                this.projectileDamage = this.originalProjectileDamage; // Revert
                System.out.println(this.name + " attack boost wore off. Projectile stats reverted.");
            }
        }
    }
    public final void step(Game game) {// DONT CHANGE
        if(!isAlive()) {
            return; // If the robot is dead, do nothing
        }
        updatePowerUpEffects();

        if(Math.abs(xMovement) + Math.abs(yMovement) > 1) {
            throw new IllegalArgumentException("You can only move in one direction at a time, use xMovement and yMovement to set the direction");
        }
        // shoot
        if (shoot && canAttack()) {
            Projectile p = new Projectile(x + Utilities.ROBOT_SIZE / 2 - Utilities.PROJECTILE_SIZE / 2, y + Utilities.ROBOT_SIZE / 2 - Utilities.PROJECTILE_SIZE / 2, xTarget, yTarget, projectileSpeed, projectileDamage, projectileImage,
                    this);
            game.addProjectile(p);
            attackCurCooldown = attackMaxCooldown;
        }

        // --- MOVEMENT ---
        int effectiveSpeed = this.speed; // this.speed is int
        if (isTileMud(this.x, this.y, game.getMap())) {
            effectiveSpeed /= 2.0;
        }

        for (int i = 0; i < effectiveSpeed; i++) { // Iterate 'effectiveSpeedSteps' times
            int potentialNextX = x + xMovement;
            int potentialNextY = y + yMovement; // Y remains unchanged for X movement
            if (canMoveTo(potentialNextX, potentialNextY, game)) {
                x = potentialNextX;
                y = potentialNextY; // Update Y only if X movement is successful
            } else {
                break; // Collision detected, stop moving in X
            }
        }

        // --- END MOVEMENT ---

        if (attackCurCooldown > 0) {
            attackCurCooldown--;
        }

        // clear out the things that should be changed in think
        xMovement = 0;
        yMovement = 0;
        xTarget = 0;
        yTarget = 0;
        shoot = false;
    }

    public int getHealth() {
        return health;
    }

    public int getMaxHealth() {
        return maxHealth;
    }

    public int getSpeed() {
        return speed;
    }

    public String getName() {
        return name;
    }

    public BufferedImage getImage() {
        return image;
    }

    public BufferedImage getProjectileImage() {
        return projectileImage;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public int getHealthPoints() {
        return healthPoints;
    }

    public int getSpeedPoints() {
        return speedPoints;
    }

    public int getAttackSpeedPoints() {
        return attackSpeedPoints;
    }

    public int getProjectileStrengthPoints() {
        return projectileStrengthPoints;
    }

    public void takeDamage(int amount) {
        this.health -= amount;
        if (this.health < 0) {
            this.health = 0;
        }
    }

    public boolean isAlive() {
        return this.health > 0;
    }

    public void setSuccessfulThink(boolean successful) {
       
        this.successfulThink = successful;
    }

    public boolean isSuccessfulThink() {
        return successfulThink;
    }

    // Getters for power-up effects
    public boolean hasSpeedBoost() {
        return speedBoostDuration > 0;
    }

    public boolean hasAttackBoost() {
        return attackBoostDuration > 0;
    }
}