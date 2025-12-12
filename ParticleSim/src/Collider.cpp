#include "Collider.h"
#include "Config.h"
#include <algorithm>
#include <cmath>

static float W_poly6(float r, float h) {
  if (r >= 0 && r <= h) {
    float hr2 = h * h - r * r;
    return 315.0f / (64.0f * M_PI * pow(h, 9)) * hr2 * hr2 * hr2;
  }
  return 0.0f;
}

static void gradW_spiky(float rx, float ry, float rz, float r, float h,
                        float &gx, float &gy, float &gz) {
  if (r > 0 && r <= h) {
    float coeff = -45.0f / (M_PI * pow(h, 6)) * (h - r) * (h - r) / r;
    gx = coeff * rx;
    gy = coeff * ry;
    gz = coeff * rz;
  } else {
    gx = gy = gz = 0.0f;
  }
}

static float laplacianW_viscosity(float r, float h) {
  if (r > 0 && r <= h) {
    return (45.0f / (M_PI * pow(h, 6))) * (h - r);
  }
  return 0.0f;
}

void Collider::isOutOfBounds(Particle &particle) {
  if (isXOutOfBounds(particle))
    reverseXVelocity(particle);
  if (isYOutOfBounds(particle))
    reverseYVelocity(particle);
  if (isZOutOfBounds(particle))
    reverseZVelocity(particle);
}

bool Collider::isXOutOfBounds(Particle &particle) {
  if (particle.x < Config::xmin) {
    particle.x = Config::xmin;
    return true;
  } else if (particle.x > Config::xmax) {
    particle.x = Config::xmax;
    return true;
  }
  return false;
}

bool Collider::isYOutOfBounds(Particle &particle) {
  if (particle.y < Config::ymin) {
    particle.y = Config::ymin;
    return true;
  } else if (particle.y > Config::ymax) {
    particle.y = Config::ymax;
    return true;
  }
  return false;
}

bool Collider::isZOutOfBounds(Particle &particle) {
  if (particle.z < Config::zmin) {
    particle.z = Config::zmin;
    return true;
  } else if (particle.z > Config::zmax) {
    particle.z = Config::zmax;
    return true;
  }
  return false;
}

void Collider::reverseXVelocity(Particle &particle) { particle.vx *= -0.5f; }
void Collider::reverseYVelocity(Particle &particle) { particle.vy *= -0.5f; }
void Collider::reverseZVelocity(Particle &particle) { particle.vz *= -0.5f; }

// --- HARDE BOTSINGEN (BACKUP) ---
// Dit lost extreme overlaps op door positie direct aan te passen
void Collider::resolveParticleCollision(Particle &A, Particle &B) {
  float dx = B.x - A.x;
  float dy = B.y - A.y;
  float dz = B.z - A.z;
  float dist2 = dx * dx + dy * dy + dz * dz;
  float minDist = Config::diameter; // Harde shell

  if (dist2 >= minDist * minDist)
    return;

  float dist = std::sqrt(dist2);
  if (dist < 0.0001f)
    dist = 0.0001f;

  float nx = dx / dist;
  float ny = dy / dist;
  float nz = dz / dist;

  float overlap = (minDist - dist) * 0.5f;

  // Positie correctie
  A.x -= nx * overlap;
  A.y -= ny * overlap;
  A.z -= nz * overlap;
  B.x += nx * overlap;
  B.y += ny * overlap;
  B.z += nz * overlap;

  // Velocity correctie (simpele bounce)
  float dvx = B.vx - A.vx;
  float dvy = B.vy - A.vy;
  float dvz = B.vz - A.vz;
  float relVel = dvx * nx + dvy * ny + dvz * nz;

  if (relVel > 0)
    return;

  float j = -(1.0f + 0.5f) * relVel * 0.5f;
  A.vx -= j * nx;
  A.vy -= j * ny;
  A.vz -= j * nz;
  B.vx += j * nx;
  B.vy += j * ny;
  B.vz += j * nz;
}

void Collider::checkCollions(const Grid &grid,
                             std::vector<Particle> &particles) {
  const auto &cells = grid.getCells();
  int Nx = grid.getNx();
  int Ny = grid.getNy();
  int Nz = grid.getNz();
  for (int cz = 0; cz < Nz; ++cz) {
    for (int cy = 0; cy < Ny; ++cy) {
      for (int cx = 0; cx < Nx; ++cx) {
        int cellID = cx + cy * Nx + cz * Nx * Ny;
        for (int dz = -1; dz <= 1; ++dz) {
          for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = cx + dx;
              int ny = cy + dy;
              int nz = cz + dz;
              if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 ||
                  nz >= Nz)
                continue;
              int neighborID = nx + ny * Nx + nz * Nx * Ny;
              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  if (i < j)
                    resolveParticleCollision(particles[i], particles[j]);
                }
              }
            }
          }
        }
      }
    }
  }
}

void Collider::applyPressure(const Grid &grid, std::vector<Particle> &particles,
                             float dt) {
  int N = particles.size();
  float h = Config::smoothingRadius;
  float mass = Config::particleMass;
  float k = Config::stiffness;
  float rho0 = Config::restDensity;
  float mu = Config::viscosity;
  float repulsionK = Config::repulsionStiffness;
  float maxF = Config::maxPressureForce;

  std::vector<float> rho(N, 0.0f);
  std::vector<float> pressure(N, 0.0f);
  std::vector<float> fx(N, 0), fy(N, 0), fz(N, 0);

  const auto &cells = grid.getCells();
  int Nx = grid.getNx(), Ny = grid.getNy(), Nz = grid.getNz();

  // Density Loop
  for (int cz = 0; cz < Nz; ++cz)
    for (int cy = 0; cy < Ny; ++cy)
      for (int cx = 0; cx < Nx; ++cx) {
        int cellID = cx + cy * Nx + cz * Nx * Ny;
        for (int dz = -1; dz <= 1; ++dz)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 ||
                  nz >= Nz)
                continue;
              int neighborID = nx + ny * Nx + nz * Nx * Ny;

              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  // Geen if (i==j) continue hier!
                  float rx = particles[i].x - particles[j].x;
                  float ry = particles[i].y - particles[j].y;
                  float rz = particles[i].z - particles[j].z;
                  float r2 = rx * rx + ry * ry + rz * rz;

                  if (r2 < h * h) {
                    rho[i] += mass * W_poly6(std::sqrt(r2), h);
                  }
                }
              }
            }
      }

  // Pressure Calculation (Non-negative)
  for (int i = 0; i < N; ++i) {
    if (rho[i] < 0.001f)
      rho[i] = 0.001f; // Veiligheid
    pressure[i] = std::max(0.0f, k * (rho[i] - rho0));
  }

  // Force Loop (Pressure + Viscosity + Repulsion)
  for (int cz = 0; cz < Nz; ++cz)
    for (int cy = 0; cy < Ny; ++cy)
      for (int cx = 0; cx < Nx; ++cx) {
        int cellID = cx + cy * Nx + cz * Nx * Ny;
        for (int dz = -1; dz <= 1; ++dz)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = cx + dx, ny = cy + dy, nz = cz + dz;
              if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 ||
                  nz >= Nz)
                continue;
              int neighborID = nx + ny * Nx + nz * Nx * Ny;

              for (int i : cells[cellID]) {
                for (int j : cells[neighborID]) {
                  if (i >= j)
                    continue; // Single check pairs

                  float rx = particles[i].x - particles[j].x;
                  float ry = particles[i].y - particles[j].y;
                  float rz = particles[i].z - particles[j].z;
                  float r2 = rx * rx + ry * ry + rz * rz;

                  // Sla over als ze te ver zijn
                  if (r2 >= h * h || r2 < 0.000001f)
                    continue;

                  float r = std::sqrt(r2);

                  // --- A. SPH KRACHTEN (Druk + Viscositeit) ---
                  float gx, gy, gz;
                  gradW_spiky(rx, ry, rz, r, h, gx, gy, gz);

                  // 1. Druk (Pressure) - CLAMPED
                  float pressTerm = -mass * mass *
                                    (pressure[i] / (rho[i] * rho[i]) +
                                     pressure[j] / (rho[j] * rho[j]));

                  if (pressTerm > maxF)
                    pressTerm = maxF;
                  if (pressTerm < -maxF)
                    pressTerm = -maxF;

                  float f_press_x = pressTerm * gx;
                  float f_press_y = pressTerm * gy;
                  float f_press_z = pressTerm * gz;

                  // 2. Viscositeit (Demping)
                  float vx_diff = particles[j].vx - particles[i].vx;
                  float vy_diff = particles[j].vy - particles[i].vy;
                  float vz_diff = particles[j].vz - particles[i].vz;

                  float laplacian = laplacianW_viscosity(r, h);
                  float viscTerm = mu * mass * (1.0f / rho[j]) * laplacian;

                  float f_visc_x = viscTerm * vx_diff;
                  float f_visc_y = viscTerm * vy_diff;
                  float f_visc_z = viscTerm * vz_diff;

                  // REPULSION KRACHT
                  // Werkt alleen als ze elkaar fysiek raken (r < diameter)
                  float f_repulse_x = 0, f_repulse_y = 0, f_repulse_z = 0;

                  if (r < Config::diameter) {
                    float penetration = Config::diameter - r;
                    float repulseForce =
                        repulsionK * penetration; // Hooke's law achtig

                    // Normaal vector
                    float nx = rx / r;
                    float ny = ry / r;
                    float nz = rz / r;

                    f_repulse_x = repulseForce * nx;
                    f_repulse_y = repulseForce * ny;
                    f_repulse_z = repulseForce * nz;
                  }

                  // TOTALE KRACHT OPTELLEN
                  float total_fx = f_press_x + f_visc_x + f_repulse_x;
                  float total_fy = f_press_y + f_visc_y + f_repulse_y;
                  float total_fz = f_press_z + f_visc_z + f_repulse_z;

                  fx[i] += total_fx;
                  fy[i] += total_fy;
                  fz[i] += total_fz;

                  fx[j] -= total_fx;
                  fy[j] -= total_fy;
                  fz[j] -= total_fz;
                }
              }
            }
      }

  // Update Velocities
  for (int i = 0; i < N; ++i) {
    particles[i].vx += dt * fx[i] / mass;
    particles[i].vy += dt * fy[i] / mass;
    particles[i].vz += dt * fz[i] / mass;
  }
}
