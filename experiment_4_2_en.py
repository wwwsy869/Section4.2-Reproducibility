import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import os
import openpyxl

def process_matrices(original_matrices, num_samples, sheet_name_prefix):
    if not original_matrices:
        return []
    num_rows = original_matrices[0].shape[0]
    num_periods = len(original_matrices)
    perturbation_matrices = [np.random.normal(0, 0, (num_rows, num_samples)) for _ in range(num_periods)]
    new_matrices = [np.zeros((num_rows, num_samples)) for _ in range(num_periods)]
    for r in range(num_rows):
        for c in range(num_samples):
            selected_indices = np.random.choice(list(range(num_periods)), 2, replace=False)
            unselected_index = list(set(range(num_periods)) - set(selected_indices))[0]
            for index in selected_indices:
                rho = perturbation_matrices[index][r, c]
                new_matrices[index][r, c] = original_matrices[index][r, c] * (1 + rho)
            rho1, rho2 = perturbation_matrices[selected_indices[0]][r, c], perturbation_matrices[selected_indices[1]][
                r, c]
            val_t1, val_t2 = original_matrices[selected_indices[0]][r, c], original_matrices[selected_indices[1]][r, c]
            new_val_unselected = original_matrices[unselected_index][r, c] - (val_t1 * rho1 + val_t2 * rho2) / 2
            new_matrices[unselected_index][r, c] = max(new_val_unselected, 0)
    try:
        print(f"  Saving sample data to: {sheet_name_prefix}.xlsx ...")
        with pd.ExcelWriter(f"{sheet_name_prefix}.xlsx", engine='openpyxl') as wb:
            original_df = pd.DataFrame([np.array(matrix).flatten() for matrix in original_matrices]).T
            original_df.columns = [f'M{i + 1}' for i in range(len(original_matrices))]
            original_df.to_excel(wb, sheet_name=f"{sheet_name_prefix}_original", index=False)
            new_df = pd.DataFrame([np.array(matrix).flatten() for matrix in new_matrices]).T
            new_df.columns = [f'M{i + 1}' for i in range(len(new_matrices))]
            new_df.to_excel(wb, sheet_name=f"{sheet_name_prefix}_uncertain", index=False)
            perturbation_df = pd.DataFrame([np.array(matrix).flatten() for matrix in perturbation_matrices]).T
            perturbation_df.columns = [f'rho_{i + 1}' for i in range(len(perturbation_matrices))]
            perturbation_df.to_excel(wb, sheet_name=f"{sheet_name_prefix}_perturb", index=False)
        print(f"  Successfully saved: {sheet_name_prefix}.xlsx")
    except Exception as e:
        print(f"Error saving {sheet_name_prefix}.xlsx: {e}")
    return new_matrices

class BaseDROModel:
    def __init__(self, name="Base_DRO_Model"):
        self.name = name
        self.I, self.J, self.T = 0, 0, 0
        self.d, self.D1, self.D2, self.c, self.c_prime = [None] * 5
        self.r, self.b, self.b_prime, self.I0, self.I0_prime = [None] * 5
        self.tc, self.tc_prime, self.dc = 0, 0, 0
        self.H, self.NV, self.NH, self.VCap_T, self.VCap_H = [None] * 5
        self.alpha = None
        self.delta, self.vartheta, self.a_bar, self.a_underline = [None] * 4
        self.mu, self.sigma, self.d_prime_bar, self.d_prime_underline = [None] * 4
        self.epsilon, self.rho, self.h_bar, self.h_underline = [None] * 4
        self.model = None

    def solve(self):
        self.model.setParam('OutputFlag', 0)
        self.model.setParam('TimeLimit', 300)
        self.model.optimize()

        if self.model.status == GRB.OPTIMAL:
            print(f"Model {self.name} solved to optimality.")
            obj_val = self.model.ObjVal

            shortage_cost = sum(
                self.r[j, t] * self.model.getVarByName(f"s[{j},{t}]").X
                for j in range(self.J) for t in range(self.T)
            )

            delay_cost_val = sum(
                self.epsilon[i, j, t] * self.model.getVarByName(f"phi1[{i},{j},{t}]").X +
                self.rho[i, j, t] * self.model.getVarByName(f"phi2[{i},{j},{t}]").X +
                self.model.getVarByName(f"phi3[{i},{j},{t}]").X
                for i in range(self.I) for j in range(self.J) for t in range(self.T)
            )

            transshipment_cost_val = 0
            for t in range(self.T):
                for i in range(self.I):
                    for i_prime in range(self.I):
                        if i == i_prime: continue
                        y_val = self.model.getVarByName(f"y[{i},{i_prime},{t}]").X
                        transshipment_cost_val += self.tc * self.D2[i, i_prime] * y_val
                        z_val = self.model.getVarByName(f"z[{i},{i_prime},{t}]").X
                        transshipment_cost_val += self.tc_prime * self.D2[i, i_prime] * z_val

            direct_transport_cost_val = obj_val - shortage_cost - delay_cost_val - transshipment_cost_val

            print(f"  Total cost: {obj_val:.2f}")
            print(f"  Direct transportation cost: {direct_transport_cost_val:.2f} (truck + helicopter)")
            print(f"  Lateral transshipment cost: {transshipment_cost_val:.2f}")
            print(f"  Delay cost: {delay_cost_val:.2f}")
            print(f"  Shortage cost: {shortage_cost:.2f}")
            print("-" * 30)

            return {
                "Total cost": obj_val,
                "Direct Transport cost": direct_transport_cost_val,
                "Transshipment cost": transshipment_cost_val,
                "Shortage cost": shortage_cost,
                "Delay cost": delay_cost_val
            }
        else:
            print(f"Model {self.name} did not reach optimality. Status: {self.model.status}")
            return None

    def _build_common_model(self, y_ub=GRB.INFINITY, z_ub=GRB.INFINITY, g_ub=GRB.INFINITY, p_ub=GRB.INFINITY, nh_ub=GRB.INFINITY):
        self.model = gp.Model(self.name)
        q = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="q")
        g = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, ub=g_ub, name="g")
        p = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, ub=p_ub, name="p")
        y = self.model.addVars(self.I, self.I, self.T, vtype=GRB.CONTINUOUS, lb=0, ub=y_ub, name="y")
        z = self.model.addVars(self.I, self.I, self.T, vtype=GRB.CONTINUOUS, lb=0, ub=z_ub, name="z")
        s = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="s")
        Inv = self.model.addVars(self.I, self.T, vtype=GRB.CONTINUOUS, lb=0, name="I")
        Inv_prime = self.model.addVars(self.I, self.T, vtype=GRB.CONTINUOUS, lb=0, name="I_prime")
        nv = self.model.addVars(self.I, self.J, self.T, vtype=GRB.INTEGER, lb=0, name="nv")
        nh = self.model.addVars(self.I, self.J, self.T, vtype=GRB.INTEGER, lb=0, ub=nh_ub, name="nh")
        eta = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="eta")
        xi = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="xi")
        varpi = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="varpi")
        omega1 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="omega1")
        omega2 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="omega2")
        omega3 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="omega3")
        omega4 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="omega4")
        phi1 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="phi1")
        phi2 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="phi2")
        phi3 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="phi3")
        beta1 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="beta1")
        beta2 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="beta2")
        beta3 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="beta3")
        beta4 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="beta4")
        nu1 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="nu1")
        nu2 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="nu2")
        nu3 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="nu3")
        nu4 = self.model.addVars(self.I, self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="nu4")
        gamma_bar0 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="gamma_bar0")
        gamma_bar1 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, name="gamma_bar1")
        gamma_bar2 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="gamma_bar2")
        lambda_bar = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="lambda_bar")
        theta_bar1 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar1")
        theta_bar2 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar2")
        theta_bar3 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar3")
        theta_bar4 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar4")
        theta_bar5 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar5")
        theta_bar6 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar6")
        theta_bar7 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar7")
        theta_bar8 = self.model.addVars(self.J, self.T, vtype=GRB.CONTINUOUS, lb=0, name="theta_bar8")
        for t in range(self.T):
            for i in range(self.I):
                I_prev = self.I0[i] if t == 0 else Inv[i, t - 1]
                I_prime_prev = self.I0_prime[i] if t == 0 else Inv_prime[i, t - 1]
                self.model.addConstr(Inv[i, t] == I_prev + self.b[i, t] + gp.quicksum(
                    y[i_prime, i, t] for i_prime in range(self.I) if i_prime != i) - gp.quicksum(
                    q[i, j, t] for j in range(self.J)) - gp.quicksum(
                    y[i, i_prime, t] for i_prime in range(self.I) if i_prime != i))
                self.model.addConstr(Inv_prime[i, t] == I_prime_prev + self.b_prime[i, t] + gp.quicksum(
                    z[i_prime, i, t] for i_prime in range(self.I) if i_prime != i) - gp.quicksum(
                    p[i, j, t] + g[i, j, t] for j in range(self.J)) - gp.quicksum(
                    z[i, i_prime, t] for i_prime in range(self.I) if i_prime != i))
                self.model.addConstr(gp.quicksum(q[i, j, t] for j in range(self.J)) + gp.quicksum(
                    y[i, i_prime, t] for i_prime in range(self.I) if i_prime != i) <= I_prev + self.b[
                                         i, t] + gp.quicksum(
                    y[i_prime, i, t] for i_prime in range(self.I) if i_prime != i))
                self.model.addConstr(gp.quicksum(p[i, j, t] + g[i, j, t] for j in range(self.J)) + gp.quicksum(
                    z[i, i_prime, t] for i_prime in range(self.I) if i_prime != i) <= I_prime_prev + self.b_prime[
                                         i, t] + gp.quicksum(
                    z[i_prime, i, t] for i_prime in range(self.I) if i_prime != i))
            for j in range(self.J):
                self.model.addConstr(gp.quicksum(q[i, j, t] for i in range(self.I)) + s[j, t] == self.d[j, t])
            self.model.addConstr(gp.quicksum(nv[i, j, t] for i in range(self.I) for j in range(self.J)) <= self.NV[t])
            self.model.addConstr(gp.quicksum(nh[i, j, t] for i in range(self.I) for j in range(self.J)) <= self.NH[t])
        for i, j, t in q.keys():
            self.model.addConstr(q[i, j, t] + g[i, j, t] <= self.VCap_T * nv[i, j, t])
            self.model.addConstr(p[i, j, t] <= self.VCap_H * nh[i, j, t])
            self.model.addConstr(
                varpi[i, j, t] >= self.a_bar[i, j, t] * omega1[i, j, t] - self.a_underline[i, j, t] * omega2[
                    i, j, t] + (omega3[i, j, t] - omega4[i, j, t]) * self.delta[i, j, t])
            self.model.addConstr(
                omega1[i, j, t] - omega2[i, j, t] + omega3[i, j, t] - omega4[i, j, t] >= self.c[i, j] * (
                            q[i, j, t] + g[i, j, t]) - eta[i, j, t])
            self.model.addConstr(omega3[i, j, t] + omega4[i, j, t] <= xi[i, j, t])
            self.model.addConstr(
                phi3[i, j, t] >= self.h_bar[i, j, t] * beta1[i, j, t] - self.h_underline[i, j, t] * beta2[i, j, t] +
                self.epsilon[i, j, t] * (beta3[i, j, t] - beta4[i, j, t]))
            self.model.addConstr(beta1[i, j, t] - beta2[i, j, t] + beta3[i, j, t] - beta4[i, j, t] >= -phi1[i, j, t])
            self.model.addConstr(beta3[i, j, t] + beta4[i, j, t] <= phi2[i, j, t])
            self.model.addConstr(
                phi3[i, j, t] >= self.h_bar[i, j, t] * nu1[i, j, t] - self.h_underline[i, j, t] * nu2[i, j, t] +
                self.epsilon[i, j, t] * (nu3[i, j, t] - nu4[i, j, t]) - self.H[j, t] * self.dc * g[i, j, t])
            self.model.addConstr(
                nu1[i, j, t] - nu2[i, j, t] + nu3[i, j, t] - nu4[i, j, t] >= self.dc * g[i, j, t] - phi1[i, j, t])
            self.model.addConstr(nu3[i, j, t] + nu4[i, j, t] <= phi2[i, j, t])
        for j, t in s.keys():
            self.model.addConstr(
                gamma_bar0[j, t] + gamma_bar1[j, t] * self.mu[j, t] - gamma_bar2[j, t] * self.sigma[j, t] - (
                            1 - self.alpha[j, t]) * lambda_bar[j, t] >= 0)
            self.model.addConstr(
                gamma_bar0[j, t] + (theta_bar1[j, t] - theta_bar2[j, t]) * self.mu[j, t] + theta_bar3[j, t] *
                self.d_prime_bar[j, t] - theta_bar4[j, t] * self.d_prime_underline[j, t] - lambda_bar[j, t] <= 0)
            self.model.addConstr(
                gamma_bar1[j, t] - theta_bar1[j, t] + theta_bar2[j, t] - theta_bar3[j, t] + theta_bar4[j, t] <= 0)
            self.model.addConstr(theta_bar1[j, t] + theta_bar2[j, t] - gamma_bar2[j, t] <= 0)
            self.model.addConstr(
                gamma_bar0[j, t] + (theta_bar5[j, t] - theta_bar6[j, t]) * self.mu[j, t] + theta_bar7[j, t] *
                self.d_prime_bar[j, t] - theta_bar8[j, t] * self.d_prime_underline[j, t] - gp.quicksum(
                    p[i, j, t] + g[i, j, t] for i in range(self.I)) <= 0)
            self.model.addConstr(
                gamma_bar1[j, t] - theta_bar5[j, t] + theta_bar6[j, t] - theta_bar7[j, t] + theta_bar8[j, t] + 1 <= 0)
            self.model.addConstr(theta_bar5[j, t] + theta_bar6[j, t] - gamma_bar2[j, t] <= 0)
        helocopter_cost = gp.quicksum(self.c_prime[i, j] * p[i, j, t] for i, j, t in p.keys())
        transshipment_cost = gp.quicksum(
            self.tc * self.D2[i, i_prime] * y[i, i_prime, t] for i, i_prime, t in y.keys() if i != i_prime) + \
                             gp.quicksum(
                                 self.tc_prime * self.D2[i, i_prime] * z[i, i_prime, t] for i, i_prime, t in z.keys() if
                                 i != i_prime)
        shortage_cost = gp.quicksum(self.r[j, t] * s[j, t] for j, t in s.keys())
        worst_case_transport_cost = gp.quicksum(
            self.delta[i, j, t] * eta[i, j, t] + self.vartheta[i, j, t] * xi[i, j, t] + varpi[i, j, t] for i, j, t in
            q.keys())
        worst_case_delay_cost = gp.quicksum(
            self.epsilon[i, j, t] * phi1[i, j, t] + self.rho[i, j, t] * phi2[i, j, t] + phi3[i, j, t] for i, j, t in
            g.keys())
        self.model.setObjective(
            helocopter_cost + transshipment_cost + shortage_cost + worst_case_transport_cost + worst_case_delay_cost,
            GRB.MINIMIZE)


class Model1(BaseDROModel):
    def build_model(self): self._build_common_model()


class Model2(BaseDROModel):
    def build_model(self): self._build_common_model(y_ub=0, z_ub=0)


class Model4(BaseDROModel):
    def build_model(self): self._build_common_model(p_ub=0, nh_ub=0)


class Model6(BaseDROModel):
    def build_model(self): self._build_common_model(y_ub=0, z_ub=0, p_ub=0, nh_ub=0)

if __name__ == "__main__":
    I, J, T, SAMPLES = 3, 10, 3, 100
    np.random.seed(42)
    params = {
        'd': np.random.uniform(200, 1000, (J, T)),
        'r': np.random.uniform(10, 50, (J, T)),
        'H': np.zeros((J, T)),
        'alpha': np.full((J, T), 0.05),
        'tc': 0.008, 'tc_prime': 0.02, 'dc': 0.2,
        'NV': [15] * T, 'NH': [6] * T,
        'VCap_T': 750, 'VCap_H': 500
    }
    params['D1'] = np.array([
        [50.6, 30.7, 84.3, 137.2, 69.6, 103, 111.9, 135.0, 103.6, 116.5],
        [60.7, 63.6, 41.4, 102.5, 146.8, 146.8, 57.8, 41.4, 49.5, 59.1],
        [111.6, 114.5, 77.4, 156, 197.7, 63.4, 96.4, 131.8, 41.2, 57.2]
    ])
    params['D2'] = np.array([
        [0, 94.1, 135.7],
        [94.1, 0, 90.6],
        [135.7, 90.6, 0]
    ])
    params['c'] = 0.02 * params['D1']
    params['c_prime'] = 0.1 * params['D1'] + 1.5
    params['b'] = np.array([
        [0, 2000, 4000],
        [50, 300, 2000],
        [500, 1500, 1500]
    ])

    params['b_prime'] = np.array([
        [0, 3000, 3000],
        [20, 20, 1000],
        [800, 1500, 1500]
    ])

    params['I0'] = np.array([5000, 2000, 2500])
    params['I0_prime'] = np.array([2500, 1500, 1000])

    for j in range(J):
        params['H'][j, :] = 2 * np.max(params['D1'][:, j])

    dro_params = {}

    all_samples_a, all_samples_h = [], []

    for i in range(I):
        original_a_for_i = [np.random.uniform(1, 10, (J, SAMPLES)) for t in range(T)]
        original_h_for_i_list = []
        for t in range(T):
            h_array_t = np.zeros((J, SAMPLES))
            for j in range(J):
                low = params['D1'][i, j]
                high = 5 * params['D1'][i, j]
                h_array_t[j, :] = np.random.uniform(low, high, SAMPLES)
            original_h_for_i_list.append(h_array_t)
        original_h_for_i = original_h_for_i_list

        samples_a_i_list = process_matrices(original_a_for_i, SAMPLES,
                                            sheet_name_prefix=f"road condition_{i + 1}")
        all_samples_a.append(np.stack(samples_a_i_list, axis=1))

        samples_h_i_list = process_matrices(original_h_for_i, SAMPLES,
                                            sheet_name_prefix=f"transportation time_{i + 1}")
        all_samples_h.append(np.stack(samples_h_i_list, axis=1))

    final_samples_a = np.stack(all_samples_a, axis=0)
    final_samples_h = np.stack(all_samples_h, axis=0)

    dro_params['delta'] = np.mean(final_samples_a, axis=3)
    dro_params['vartheta'] = np.mean(np.abs(final_samples_a - dro_params['delta'][:, :, :, np.newaxis]), axis=3)
    dro_params['a_bar'] = np.max(final_samples_a, axis=3)
    dro_params['a_underline'] = np.min(final_samples_a, axis=3)

    dro_params['epsilon'] = np.mean(final_samples_h, axis=3)
    dro_params['rho'] = np.mean(np.abs(final_samples_h - dro_params['epsilon'][:, :, :, np.newaxis]), axis=3)
    dro_params['h_bar'] = np.max(final_samples_h, axis=3)
    dro_params['h_underline'] = np.min(final_samples_h, axis=3)

    original_d_prime = [np.random.uniform(100, 500, (J, SAMPLES)) for t in range(T)]

    samples_d_prime_list = process_matrices(original_d_prime, SAMPLES,
                                            sheet_name_prefix="emergency demand")

    samples_d_prime = np.stack(samples_d_prime_list, axis=1)

    dro_params['mu'] = np.mean(samples_d_prime, axis=2)
    dro_params['sigma'] = np.mean(np.abs(samples_d_prime - dro_params['mu'][:, :, np.newaxis]), axis=2)
    dro_params['d_prime_bar'] = np.max(samples_d_prime, axis=2)
    dro_params['d_prime_underline'] = np.min(samples_d_prime, axis=2)
    models_to_run = {
        "LT-TMM": Model1(name="Model1_LT-TMM"),
        "TMM": Model2(name="Model2_TMM"),
        "LT-ERT-T": Model4(name="Model4_LT-ERT-T"),
        "NSU-ERT-T": Model6(name="Model6_NSU-ERT-T"),
    }

    results = {}

    for model_key, model_instance in models_to_run.items():
        print(f"\n>>> Solving model: {model_key} ({model_instance.name})")
        model_instance.I, model_instance.J, model_instance.T = I, J, T
        for p_name, p_value in params.items(): setattr(model_instance, p_name, p_value)
        for p_name, p_value in dro_params.items(): setattr(model_instance, p_name, p_value)
        model_instance.build_model()
        result = model_instance.solve()
        if result: results[model_key] = result

    if results:
        results_df = pd.DataFrame.from_dict(results, orient='index')
        print(results_df.to_string())

        def print_benefit_analysis(our_res, bench_res, benchmark_desc):
            if not bench_res: return

            investment = our_res["Transshipment cost"] - bench_res.get("Transshipment cost", 0)
            delay_reduction = bench_res["Delay cost"] - our_res["Delay cost"]
            shortage_reduction = bench_res["Shortage cost"] - our_res["Shortage cost"]
            effectiveness_gain = delay_reduction + shortage_reduction
            efficiency_gain = bench_res["Direct Transport cost"] - our_res["Direct Transport cost"]
            net_benefit = effectiveness_gain + efficiency_gain - investment

            print(f"\nComparison: Proposed model (LT-TMM) vs {benchmark_desc}")
            print("-" * 80)
            print(f"1. Transshipment investment: {investment:+.2f}")
            print(f"2. Effectiveness gain (shortage and delay reduction): {effectiveness_gain:+.2f}")
            print(f"3. Efficiency gain (direct transportation): {efficiency_gain:+.2f}")
            print(f"Net spatial-capacity synergy benefit: {net_benefit:+.2f}")
            print("=" * 80)

        if "LT-TMM" in results:
            our = results["LT-TMM"]
            print_benefit_analysis(our, results.get("TMM"), "TMM")
            print_benefit_analysis(our, results.get("LT-ERT-T"), "LT-ERT-T")
            print_benefit_analysis(our, results.get("NSU-ERT-T"), "NSU-ERT-T")
        else:
            print("Unable to calculate the spatial-capacity synergy benefit.")
        print(results_df.to_string())

        output_filename = "experiment_1_results.xlsx"
        print(f"\nSaving results and parameters to: {output_filename} ...")

        try:
            with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:

                results_df.to_excel(writer, sheet_name="Experiment_1_Results")
                scalar_params = {}
                array_params_to_save = {}

                for key, value in params.items():
                    if isinstance(value, np.ndarray):
                        array_params_to_save[f"Param_{key}"] = value
                    else:
                        scalar_params[key] = str(value)

                pd.Series(scalar_params, name="Value").to_excel(writer, sheet_name="Parameters_Scalar")

                for sheet_name, arr in array_params_to_save.items():
                    if arr.ndim == 3:
                        df = pd.DataFrame(arr.reshape(-1, arr.shape[-1]),
                                          columns=[f"T{t}" for t in range(arr.shape[-1])])
                    elif arr.ndim == 1:
                        df = pd.DataFrame(arr, columns=["Value"])
                    else:
                        df = pd.DataFrame(arr)
                    df.to_excel(writer, sheet_name=sheet_name)

                for key, arr in dro_params.items():
                    sheet_name = f"DRO_{key}"
                    if arr.ndim == 3:
                        df = pd.DataFrame(arr.reshape(-1, arr.shape[-1]),
                                          columns=[f"T{t}" for t in range(arr.shape[-1])])
                    else:
                        df = pd.DataFrame(arr)
                    df.to_excel(writer, sheet_name=sheet_name)

            print(f"All results and parameters were successfully saved to: {output_filename}")

        except Exception as e:
            print(f"Error saving results to Excel: {e}")