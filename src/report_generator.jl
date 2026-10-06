function generate_clinical_report(sample_id::String, result::MigrationResult, tropism_score::Float64, output_path::String)
    html_content = """
    <html>
    <head>
        <title>migrant.jl Clinical Tropism Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .high-risk { color: red; font-weight: bold; }
            .low-risk { color: green; font-weight: bold; }
        </style>
    </head>
    <body>
        <h1>migrant.jl Clinical Tropism Report: Patient $(sample_id)</h1>
        
        <h2>Genomic Migration Proof</h2>
        <p>Containment Score: $(round(result.containment_score * 100, digits=2))%</p>
        
        <h2>Virulence Profile</h2>
        <ul>
    """
    
    for (coords, coverage) in result.virulence_profile
        html_content *= "        <li>Gene $(coords[1])-$(coords[2]): $(round(coverage * 100, digits=2))% coverage</li>\n"
    end
    
    html_content *= """
        </ul>
        
        <h2>Tumor Tropism Risk</h2>
    """
    
    risk_percent = round(tropism_score * 100, digits=2)
    if tropism_score > 0.50
        html_content *= "    <p>Risk Score: <span class=\"high-risk\">$(risk_percent)% (High Risk)</span></p>\n"
    else
        html_content *= "    <p>Risk Score: <span class=\"low-risk\">$(risk_percent)% (Low Risk)</span></p>\n"
    end
    
    html_content *= """
    </body>
    </html>
    """
    
    open(io -> write(io, html_content), output_path, "w")
    return true
end
