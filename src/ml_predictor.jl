using MLJ
using DataFrames
import DecisionTree

const RFC = @load RandomForestClassifier pkg=DecisionTree verbosity=0

function train_tropism_model(df::DataFrame)
    y = coerce(df.class_label, Multiclass)
    X = select(df, Not([:sample_id, :class_label]))
    
    # Coerce features to Continuous and fill missing values with 0.0
    for col in names(X)
        if eltype(X[!, col]) >: Missing
            X[!, col] = coalesce.(X[!, col], 0.0)
        end
    end
    X = coerce(X, Count => Continuous, Continuous => Continuous)
    
    model = RFC(n_trees=100)
    mach = machine(model, X, y)
    fit!(mach, verbosity=0)
    return mach
end

function predict_tropism(mach::Machine, new_patient_df::DataFrame)
    cols_to_drop = intersect([:sample_id, :class_label], propertynames(new_patient_df))
    X_new = select(new_patient_df, Not(cols_to_drop))
    
    for col in names(X_new)
        if eltype(X_new[!, col]) >: Missing
            X_new[!, col] = coalesce.(X_new[!, col], 0.0)
        end
    end
    X_new = coerce(X_new, Count => Continuous, Continuous => Continuous)
    
    preds = predict(mach, X_new)
    prob = pdf(preds[1], 1)
    return convert(Float64, prob)
end
