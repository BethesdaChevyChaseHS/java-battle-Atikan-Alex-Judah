package bcc.javaJostle;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class MyRobot extends Robot{
    private int lastX;
    private int lastY;
    private int stuckCounter = 0;
    private int lastBotX = -1;
    private int lastBotY = -1;

    public MyRobot(int x, int y){
        super(x, y, 2, 3, 2, 3,"bob", "myRobotImage.png", "defaultProjectileImage.png");
        // Health: 3, Speed: 3, Attack Speed: 2, Projectile Strength: 2
        // Total = 10
        // Image name is "myRobotImage.png"
    }

    private class Tile {
        private int x;
        private int y;

        public Tile(int x, int y) {
            this.x = x;
            this.y = y;
        }

        public int tileX() {
            return x;
        }

        public int tileY() {
            return y;
        }
    }

    public void think(ArrayList<Robot> robots, ArrayList<Projectile> projectiles, Map map, ArrayList<PowerUp> powerups) {
        if(powerups.size() > 0) {
            dijkstras(map, projectiles, powerups.get(0).getX(), powerups.get(0).getY());
        } else {
            Robot nearest = null;
            double minDist = Double.MAX_VALUE;
            for (Robot robot : robots) {
                if (robot != this && robot.isAlive()) {
                    double dx = getX() - robot.getX();
                    double dy = getY() - robot.getY();
                    double dist = dx * dx + dy * dy;
                    if (dist < minDist) {
                        minDist = dist;
                        nearest = robot;
                    }
                }
            }
            if (nearest != null) {
                int[][] tiles = map.getTiles();
                int width = tiles[0].length;
                int height = tiles.length;
                double enemyX = nearest.getX();
                double enemyY = nearest.getY();

                double maxDist = -1;
                int farthestX = 0, farthestY = 0;
                for (int x = 0; x < width; x++) {
                    for (int y = 0; y < height; y++) {
                        if (tiles[y][x] != Utilities.WALL) {
                            double dist = Math.pow(enemyX - x * Utilities.TILE_SIZE, 2) + Math.pow(enemyY - y * Utilities.TILE_SIZE, 2);
                            if (dist > maxDist) {
                                maxDist = dist;
                                farthestX = x * Utilities.TILE_SIZE;
                                farthestY = y * Utilities.TILE_SIZE;
                            }
                        }
                    }
                }
                dijkstras(map, projectiles, farthestX, farthestY);
            }
        }
        if(canAttack()){
            for(Robot robot : robots) {
                if (robot != this && robot.isAlive() ){
                    predictiveAimAndShoot(robot);
                    break;
                }
            }
        }
        int currentX = getX();
        int currentY = getY();
        if (currentX == lastBotX && currentY == lastBotY) {
            stuckCounter++;
        } else {
            stuckCounter = 0;
        }
        lastBotX = currentX;
        lastBotY = currentY;
        if (stuckCounter >= 5) {
            int[][] tiles = map.getTiles();
            int width = tiles[0].length;
            int height = tiles.length;
            for (int tries = 0; tries < 10; tries++) {
                int rx = (int)(Math.random() * width);
                int ry = (int)(Math.random() * height);
                if (tiles[ry][rx] != Utilities.WALL) {
                    dijkstras(map, projectiles, rx * Utilities.TILE_SIZE, ry * Utilities.TILE_SIZE);
                    break;
                }
            }
            stuckCounter = 0;
        }
    }

    private void predictiveAimAndShoot(Robot target) {
        double shooterX = getX();
        double shooterY = getY();
        double targetX = target.getX();
        double targetY = target.getY();
        double targetVelX = target.xMovement * target.getSpeed();
        double targetVelY = target.yMovement * target.getSpeed();
        double projectileSpeed = getProjectileStrengthPoints() + 5;
        double relX = targetX - shooterX;
        double relY = targetY - shooterY;
        double a = targetVelX * targetVelX + targetVelY * targetVelY - projectileSpeed * projectileSpeed;
        double b = 2 * (relX * targetVelX + relY * targetVelY);
        double c = relX * relX + relY * relY;
        double discriminant = b * b - 4 * a * c;
        double interceptTime;
        if (a == 0) {
            interceptTime = -c / b;
        } else if (discriminant >= 0) {
            double sqrtDisc = Math.sqrt(discriminant);
            double t1 = (-b + sqrtDisc) / (2 * a);
            double t2 = (-b - sqrtDisc) / (2 * a);
            interceptTime = Math.min(t1, t2);
            if (interceptTime < 0) interceptTime = Math.max(t1, t2);
        } else {
            interceptTime = 0;
        }
        if (interceptTime < 0) interceptTime = 0;
        double aimX = targetX + targetVelX * interceptTime;
        double aimY = targetY + targetVelY * interceptTime;
        shootAtLocation((int) aimX, (int) aimY);
    }

    public void dijkstras(Map map, ArrayList<Projectile> projectiles, double targetX, double targetY) {
        int[][] tiles = map.getTiles();
        int width = map.getTiles()[0].length;
        int height = map.getTiles().length;
        boolean[][] visited = new boolean[width][height];
        int[][] distances = new int[width][height];
        int robotTileX = getX() / Utilities.TILE_SIZE;
        int robotTileY = getY() / Utilities.TILE_SIZE;
        int robotTileXPlus = (getX() + Utilities.ROBOT_SIZE) / Utilities.TILE_SIZE;
        int robotTileYPlus = (getY() + Utilities.ROBOT_SIZE) / Utilities.TILE_SIZE;
        int targetTileX = (int) targetX/ Utilities.TILE_SIZE;
        int targetTileY = (int) targetY / Utilities.TILE_SIZE;

        boolean[][] dangerMap = new boolean[width][height];
        for (Projectile p : projectiles) {
            if(p.getOwner() != this) {
                double px = p.getX();
                double py = p.getY();
                double angle = Math.toRadians(p.getAngle());
                double dx = Math.cos(angle);
                double dy = Math.sin(angle);

            // Predict the next 5 tiles
                for (int t = 0; t < 5; t++) {
                    int tx = (int)(px / Utilities.TILE_SIZE);
                    int ty = (int)(py / Utilities.TILE_SIZE);
                    if (tx >= 0 && tx < width && ty >= 0 && ty < height) {
                        dangerMap[tx][ty] = true;
                    }
                    px += dx * Utilities.TILE_SIZE;
                    py += dy * Utilities.TILE_SIZE;
                }
            }
        }

        for(int i = 0; i < width; i++) {
            for(int j = 0; j < height; j++) {
                if(tiles[j][i] == Utilities.WALL || dangerMap[i][j]) {  
                    visited[i][j] = true;
                }else{
                    visited[i][j] = false;
                }
                if(lastX == -1 || lastY == -1) {
                    if(robotTileXPlus == i && robotTileYPlus == j) {
                        distances[i][j] = 0;
                    }else{
                        distances[i][j] = Integer.MAX_VALUE;
                    }
                }else {
                    if(robotTileX == i && robotTileY == j) {
                        distances[i][j] = 0;
                    }else{
                        distances[i][j] = Integer.MAX_VALUE;
                    }
                }
            }
        }

        if(tiles[targetTileY][targetTileX] == Utilities.WALL || dangerMap[targetTileX][targetTileY]) {
            System.out.println(targetTileX + ", " + targetTileY);
            return;
        }

        Queue<Tile> priority = new LinkedList<>();
        if(lastX == -1 || lastY == -1) {
            priority.add(new Tile(robotTileXPlus, robotTileYPlus));
            visited[robotTileXPlus][robotTileYPlus] = true;
        } else {
            priority.add(new Tile(robotTileX, robotTileY));
            visited[robotTileX][robotTileY] = true;
        }
        while(!priority.isEmpty()) {
            int curr = priority.size();
            for(int i = 0; i < curr; i++) {
                Tile tile = priority.poll();
                int tileX = tile.tileX();
                int tileY = tile.tileY();
                if(tileX + 1 < visited.length && !visited[tileX + 1][tileY]) {
                    visited[tileX + 1][tileY] = true;
                    distances[tileX + 1][tileY] = distances[tileX][tileY] + 1;
                    priority.add(new Tile(tileX + 1, tileY));
                }
                if(tileX - 1 >= 0 &&!visited[tileX - 1][tileY]) {
                    visited[tileX - 1][tileY] = true;
                    distances[tileX - 1][tileY] = distances[tileX][tileY] + 1;
                    priority.add(new Tile(tileX - 1, tileY));
                }
                if(tileY + 1 < visited[0].length && !visited[tileX][tileY + 1]) {
                    visited[tileX][tileY + 1] = true;
                    distances[tileX][tileY + 1] = distances[tileX][tileY] + 1;
                    priority.add(new Tile(tileX, tileY + 1));
                }
                if(tileY - 1 >= 0 && !visited[tileX][tileY - 1]) {
                    visited[tileX][tileY - 1] = true;
                    distances[tileX][tileY - 1] = distances[tileX][tileY] + 1;
                    priority.add(new Tile(tileX, tileY - 1));
                }
            }
        }

        int x = (int) targetTileX;
        int y = (int) targetTileY;
        int currDist = distances[x][y];
        while(currDist != 1) {
            System.out.println(currDist + ", " + x + ", " + y);
            int min = Math.min(Math.min(x + 1 < distances.length ? distances[x + 1][y] : Integer.MAX_VALUE, x - 1 >= 0 ? distances[x - 1][y] : Integer.MAX_VALUE), Math.min(y + 1 < distances[0].length ? distances[x][y + 1] : Integer.MAX_VALUE, y - 1 >= 0 ? distances[x][y - 1] : Integer.MAX_VALUE));
            if(x + 1 < distances.length && distances[x + 1][y] == min) {
                x += 1;
            }else if(x - 1 >= 0 && distances[x - 1][y] == min) {
                x -= 1;
            }else if(y + 1 < distances[0].length && distances[x][y + 1] == min) {
                y += 1;
            }else if(y - 1 >= 0 && distances[x][y - 1] == min) {
                y -= 1;
            }
            currDist = distances[x][y];
        }
        if(lastX == -1 || lastY == -1) {
            xMovement = x - robotTileXPlus;
            yMovement = y - robotTileYPlus;
            lastX = xMovement;
            lastY = yMovement;
        } else {
            xMovement = x - robotTileX;
            yMovement = y - robotTileY;
            lastX = x - robotTileX;
            lastY = y - robotTileY;
        }
        return;
    }
}